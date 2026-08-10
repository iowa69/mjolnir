"""The only place a language model is spoken to, and the only place it may fail.

Two properties matter more than anything else this module does.

**It must not be able to stop a run.** The verdicts in a Mjolnir report are
computed in Python from stated thresholds before a model is ever contacted; the
model writes the reading over finished checks. So a host that is down, slow,
serving no model, or answering with something unparseable is a *degradation of
the prose*, never a failure of the analysis. Every method here returns an
:class:`LLMResult` carrying ``ok=False`` and a message; nothing raises into the
pipeline. That is deliberate and is the one place in Mjolnir where an error is
returned rather than raised — the caller turns it into
``Interpretation(rule_only=True, discarded_reason=...)`` and the report says so
on the page.

**It must speak whatever the operator already runs.** ollama is what sits on a
workstation; vLLM, SGLang, llama.cpp's server and TGI are what sit in a cluster,
and those expose ``POST /v1/chat/completions``. Both protocols are implemented,
the host is probed once to decide which it offers, and the answer is harvested
from whichever field the server chose to fill. That last part is not defensive
programming for its own sake: reasoning models scatter their output across
``message.content``, ``message.thinking``, ``response``, ``reasoning`` and
``reasoning_content`` depending on the server *and* the model, and reading only
the documented field is how a fleet of gates silently takes its default for a
day while every dashboard looks healthy (the lesson is tesseract-ai's, paid for
once already).

The host comes from ``$MJOLNIR_LLM_HOST``. There is a loopback default so a
laptop works with no configuration, but nothing in here assumes loopback: inside
a container the model server is a different service on a different name, the env
var is how it is named, and the port is never guessed from the backend.

Transport is :mod:`urllib.request` from the standard library rather than
``httpx``, because ``pyproject.toml`` does not carry an HTTP client and an
optional prose layer must not add a hard dependency to a resistance-calling
tool.
"""

from __future__ import annotations

import json
import os
import socket
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from mjolnir.utils import LOG

# ---------------------------------------------------------------------------
# Transport settings
#
# These are transport parameters, not scientific thresholds, so they are not in
# the config.py registry — nothing in a report is derived from them. Each still
# names why it holds the value it does.
# ---------------------------------------------------------------------------

#: SOURCE: ollama's documented default port (https://github.com/ollama/ollama,
#: `OLLAMA_HOST` defaults to 127.0.0.1:11434). Used only when neither
#: ``MJOLNIR_LLM_HOST`` nor ``OLLAMA_HOST`` is set. The literal 127.0.0.1 rather
#: than "localhost" avoids the IPv6-first resolution stall on hosts where ::1 is
#: not listened on; it is a default, and a container overrides it by env var.
DEFAULT_HOST = "http://127.0.0.1:11434"

#: Env vars consulted, in order. ``OLLAMA_HOST`` is honoured because anyone
#: running ollama in docker-compose has already set it, and making them set a
#: second variable to say the same thing invites the two to disagree.
HOST_ENV_VARS = ("MJOLNIR_LLM_HOST", "OLLAMA_HOST")

#: Bearer token for a served endpoint that wants one (vLLM's ``--api-key``,
#: llama.cpp's ``--api-key``). Absent by default; a local ollama needs none.
API_KEY_ENV_VAR = "MJOLNIR_LLM_API_KEY"

#: SOURCE: transport policy. A probe is one HTTP GET against a service that is
#: either local or on the lab network; 5 s separates "not running" from "busy"
#: without making `mjolnir doctor` feel hung.
PROBE_TIMEOUT_SECONDS = 5.0

#: SOURCE: transport policy. A 30B-class model writing a two-paragraph reading
#: on a loaded workstation GPU is a tens-of-seconds operation, and the first
#: request also pays the model load. Long enough not to abandon a working host,
#: finite because the pipeline must finish either way.
COMPLETION_TIMEOUT_SECONDS = 180.0

#: SOURCE: transport policy. The observation is capped at a few tens of KB by
#: agent/observation.py, so 32k tokens is headroom rather than a limit, and it
#: keeps the prompt cache friendly on a single-GPU host.
DEFAULT_NUM_CTX = 32768

#: Deterministic prose. The reading is a statement of fact over finished checks;
#: sampling variety in it buys nothing and makes two runs of the same sample
#: disagree in wording.
DEFAULT_TEMPERATURE = 0.0

BACKEND_OLLAMA = "ollama"
BACKEND_OPENAI = "openai"
BACKENDS = (BACKEND_OLLAMA, BACKEND_OPENAI)

#: Where each backend advertises itself, and where it takes a chat.
_PROBE_PATHS = ((BACKEND_OLLAMA, "/api/tags"), (BACKEND_OPENAI, "/v1/models"))


def host_from_env(explicit: str = "") -> str:
    """The model host: an explicit setting, then the env vars, then loopback.

    Normalising here rather than at each call site means ``ollama.lab:11434``,
    ``http://ollama.lab:11434`` and ``http://ollama.lab:8000/v1/`` all name the
    same server, which is what people actually type.
    """
    host = (explicit or "").strip()
    if not host:
        for name in HOST_ENV_VARS:
            value = os.environ.get(name, "").strip()
            if value:
                host = value
                break
    if not host:
        host = DEFAULT_HOST
    if "://" not in host:
        host = "http://" + host
    return host.rstrip("/")


def _is_loopback(host: str) -> bool:
    name = host.split("://", 1)[-1].split("/", 1)[0].split(":")[0]
    return name in ("127.0.0.1", "localhost", "::1", "[::1]", "0.0.0.0")


@dataclass
class LLMResult:
    """One exchange with the model, successful or not.

    ``ok`` False with a populated ``error`` is a first-class outcome, not an
    exception in disguise: the caller records it as the reason the report is
    rule-only. ``text`` is whatever the server actually returned, unparsed,
    because the discipline rules run on the raw answer.
    """

    text: str = ""
    ok: bool = False
    error: str = ""
    host: str = ""
    backend: str = ""
    model: str = ""
    elapsed_seconds: float = 0.0
    #: Raw attempts, oldest first, for the trace and for `mjolnir doctor -v`.
    attempts: List[str] = field(default_factory=list)

    def as_json(self) -> Optional[Dict[str, Any]]:
        """The answer parsed as a JSON object, or None if it is not one.

        None rather than ``{}``: a model that returned prose where JSON was
        asked for has not returned an empty object, and the two must stay
        distinguishable to the caller that decides whether to retry.
        """
        try:
            parsed = json.loads(self.text)
        except (ValueError, TypeError):
            return None
        return parsed if isinstance(parsed, dict) else None


def _post_json(url: str, payload: Dict[str, Any], timeout: float,
               headers: Dict[str, str]) -> Dict[str, Any]:
    """POST a JSON body and return the decoded reply.

    Raises :class:`urllib.error.URLError`, :class:`OSError` or ``ValueError``;
    every caller in this module converts those into an ``LLMResult``. Nothing
    outside this module sees them.
    """
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(url, data=body, method="POST")
    request.add_header("Content-Type", "application/json")
    for key, value in headers.items():
        request.add_header(key, value)
    opener = _opener_for(url)
    with opener.open(request, timeout=timeout) as response:
        raw = response.read().decode("utf-8", "replace")
    return json.loads(raw)


def _opener_for(url: str) -> urllib.request.OpenerDirector:
    """An opener that does not send a loopback request through a proxy.

    A workstation with ``HTTP_PROXY`` set for the outside world would otherwise
    hand ``http://127.0.0.1:11434`` to the proxy, which refuses it, and the
    report would say the model host was unreachable while ollama was running on
    the same machine. Remote hosts keep the environment's proxy settings, since
    a lab GPU server behind one is a real deployment.
    """
    if _is_loopback(url):
        return urllib.request.build_opener(urllib.request.ProxyHandler({}))
    return urllib.request.build_opener()


def _describe_transport_error(exc: BaseException) -> str:
    """A message that names what to check, not just what failed."""
    if isinstance(exc, urllib.error.HTTPError):
        detail = ""
        try:
            detail = exc.read().decode("utf-8", "replace")[:300]
        except (OSError, ValueError):
            detail = ""
        return "HTTP {0} from the model host{1}".format(
            exc.code, ": " + detail.strip() if detail.strip() else "")
    if isinstance(exc, urllib.error.URLError):
        return "cannot reach the model host: {0}".format(exc.reason)
    if isinstance(exc, socket.timeout):
        return "the model host did not answer in time"
    if isinstance(exc, ValueError):
        return "the model host returned a body that is not JSON: {0}".format(exc)
    return "{0}: {1}".format(type(exc).__name__, exc)


class LLMClient:
    """A local model, spoken to over whichever protocol it happens to offer.

    Both protocols are kept rather than collapsing onto the OpenAI one, because
    suppressing a reasoning model's thinking is worth an order of magnitude in
    latency and *every server spells the suppression differently*: ollama's
    native ``think: false`` works and its own ``/v1`` shim ignores
    ``chat_template_kwargs`` while honouring ``reasoning_effort``, and vLLM and
    SGLang are the other way round. So the OpenAI path sends both spellings, the
    native path is preferred where it exists, and :meth:`_harvest` is the
    backstop for a server that honours neither.

    The model tag is not guessed. If none is configured, the first model the
    host advertises is used and recorded in the report, because a hard-coded tag
    that was never pulled is a guaranteed 404 and "the model is unavailable" is
    then a lie about the host.
    """

    def __init__(self, model: str = "", host: str = "",
                 num_ctx: int = DEFAULT_NUM_CTX,
                 timeout: float = COMPLETION_TIMEOUT_SECONDS,
                 backend: str = "auto") -> None:
        self.host = host_from_env(host)
        self.model = (model or "").strip()
        self.num_ctx = int(num_ctx)
        self.timeout = float(timeout)
        self.requested_backend = backend
        self._backend: Optional[str] = backend if backend in BACKENDS else None
        # A host given as ".../v1" is naming the OpenAI surface explicitly.
        if self.host.endswith("/v1"):
            self.host = self.host[: -len("/v1")]
            self._backend = self._backend or BACKEND_OPENAI
        self._probe_error = ""
        self._headers: Dict[str, str] = {}
        key = os.environ.get(API_KEY_ENV_VAR, "").strip()
        if key:
            self._headers["Authorization"] = "Bearer " + key

    # -- discovery ----------------------------------------------------------

    def _get_json(self, path: str, timeout: float) -> Optional[Dict[str, Any]]:
        url = self.host + path
        request = urllib.request.Request(url, method="GET")
        for name, value in self._headers.items():
            request.add_header(name, value)
        try:
            with _opener_for(url).open(request, timeout=timeout) as response:
                if response.status != 200:
                    return None
                return json.loads(response.read().decode("utf-8", "replace"))
        except (urllib.error.URLError, OSError, ValueError) as exc:
            self._probe_error = _describe_transport_error(exc)
            return None

    def detect_backend(self) -> Optional[str]:
        """Which protocol this host speaks, or None when it is not answering.

        ollama is probed first because it answers both probes and only its
        native endpoint can switch a reasoning model's thinking off.
        """
        if self._backend:
            return self._backend
        for name, path in _PROBE_PATHS:
            if self._get_json(path, PROBE_TIMEOUT_SECONDS) is not None:
                self._backend = name
                LOG.debug("model host %s speaks %s", self.host, name)
                return name
        return None

    def available(self) -> bool:
        return self.detect_backend() is not None

    def list_models(self) -> List[str]:
        """Model tags the host is serving, in the order it lists them."""
        backend = self.detect_backend()
        if backend == BACKEND_OLLAMA:
            body = self._get_json("/api/tags", PROBE_TIMEOUT_SECONDS) or {}
            return [str(m.get("name") or m.get("model") or "")
                    for m in (body.get("models") or []) if isinstance(m, dict)]
        if backend == BACKEND_OPENAI:
            body = self._get_json("/v1/models", PROBE_TIMEOUT_SECONDS) or {}
            return [str(m.get("id") or "") for m in (body.get("data") or [])
                    if isinstance(m, dict)]
        return []

    def resolve_model(self) -> str:
        """The model tag to send, taking the host's first if none was configured."""
        if self.model:
            return self.model
        served = [m for m in self.list_models() if m]
        if served:
            self.model = served[0]
        return self.model

    def describe(self) -> str:
        """One line for `mjolnir doctor`, saying what is and is not there."""
        backend = self.detect_backend()
        if backend is None:
            return "no LLM host at {0} ({1}); the report will be rule-only".format(
                self.host, self._probe_error or "no response")
        model = self.resolve_model() or "no model served"
        return "{0} host at {1}, model {2}".format(backend, self.host, model)

    # -- completion ---------------------------------------------------------

    @staticmethod
    def _harvest(body: Dict[str, Any]) -> str:
        """The answer, from wherever this server decided to put it.

        ollama's ``/api/chat`` fills ``message.content`` and spills a reasoning
        model's output into ``message.thinking``; ``/api/generate`` uses
        ``response``; the OpenAI shape uses ``choices[0].message.content`` and
        spills into ``reasoning`` or ``reasoning_content`` depending on whose
        server it is. ``content`` is preferred everywhere it is non-empty, and
        the thinking fields are read only when it is not — an answer buried in
        the thinking channel is still an answer, and the discipline layer will
        judge it on its merits either way.
        """
        containers: List[Dict[str, Any]] = [body]
        message = body.get("message")
        if isinstance(message, dict):
            containers.append(message)
        choices = body.get("choices")
        if isinstance(choices, list) and choices and isinstance(choices[0], dict):
            containers.append(choices[0])
            inner = choices[0].get("message")
            if isinstance(inner, dict):
                containers.append(inner)
            if isinstance(choices[0].get("text"), str):
                containers.append({"content": choices[0]["text"]})
        for key in ("content", "response", "reasoning_content", "reasoning",
                    "thinking"):
            for container in containers:
                value = container.get(key)
                if isinstance(value, str) and value.strip():
                    return value.strip()
        return ""

    def _payload(self, backend: str, messages: Sequence[Dict[str, str]],
                 json_mode: bool, temperature: float) -> Dict[str, Any]:
        model = self.resolve_model()
        if backend == BACKEND_OPENAI:
            payload: Dict[str, Any] = {
                "model": model,
                "messages": list(messages),
                "temperature": temperature,
                "stream": False,
                # Both spellings of "do not think", because no server honours
                # both and each ignores the other's: ollama's /v1 shim answers
                # to `reasoning_effort`, vLLM and SGLang to
                # `chat_template_kwargs`. `_harvest` covers a server that
                # honours neither.
                "chat_template_kwargs": {"enable_thinking": False},
                "reasoning_effort": "none",
            }
            if json_mode:
                payload["response_format"] = {"type": "json_object"}
            return payload
        payload = {
            "model": model,
            "messages": list(messages),
            "stream": False,
            # Load-bearing on ollama: a reasoning model's output is routed into
            # `thinking` and `content` comes back empty, so without this the
            # answer looks like an empty reply from a healthy server.
            "think": False,
            "options": {"num_ctx": self.num_ctx, "temperature": temperature},
        }
        if json_mode:
            payload["format"] = "json"
        return payload

    def chat(self, messages: Sequence[Dict[str, str]], json_mode: bool = True,
             temperature: float = DEFAULT_TEMPERATURE) -> LLMResult:
        """One exchange. Returns a result; never raises into the pipeline."""
        started = time.time()
        backend = self.detect_backend()
        if backend is None:
            return LLMResult(
                ok=False, host=self.host, model=self.model,
                error="no LLM host answered at {0} ({1}); set {2} to point at "
                      "one, or accept a rule-only report".format(
                          self.host, self._probe_error or "no response",
                          HOST_ENV_VARS[0]),
                elapsed_seconds=time.time() - started,
            )
        if not self.resolve_model():
            return LLMResult(
                ok=False, host=self.host, backend=backend,
                error="the {0} host at {1} serves no model; pull or load one "
                      "(e.g. `ollama pull qwen3:8b`) or set {2}".format(
                          backend, self.host, "MJOLNIR_LLM_MODEL"),
                elapsed_seconds=time.time() - started,
            )

        url = self.host + ("/v1/chat/completions" if backend == BACKEND_OPENAI
                           else "/api/chat")
        payload = self._payload(backend, messages, json_mode, temperature)
        try:
            body = _post_json(url, payload, self.timeout, self._headers)
        except (urllib.error.URLError, OSError, ValueError) as exc:
            return LLMResult(ok=False, host=self.host, backend=backend,
                             model=self.model,
                             error=_describe_transport_error(exc),
                             elapsed_seconds=time.time() - started)
        text = self._harvest(body)
        elapsed = time.time() - started
        if not text:
            return LLMResult(
                ok=False, host=self.host, backend=backend, model=self.model,
                error="the model returned an empty answer",
                elapsed_seconds=elapsed, attempts=[json.dumps(body)[:300]])
        return LLMResult(text=text, ok=True, host=self.host, backend=backend,
                         model=self.model, elapsed_seconds=elapsed,
                         attempts=[text[:300]])

    def complete(self, prompt: str, system: str = "", json_mode: bool = True,
                 temperature: float = DEFAULT_TEMPERATURE) -> LLMResult:
        """A single user turn with an optional system prompt."""
        messages: List[Dict[str, str]] = []
        if system:
            messages.append({"role": "system", "content": system})
        messages.append({"role": "user", "content": prompt})
        return self.chat(messages, json_mode=json_mode, temperature=temperature)


class ScriptedClient(LLMClient):
    """An :class:`LLMClient` that answers from a list instead of a network.

    This exists so the discipline rules can be tested against exact answers —
    an invented number, a concessive clause, "susceptible" where the rule said
    no determinant was detected — without a model, a GPU or a network. It is
    also what ``--dry-run`` uses to exercise the prose path.

    It is not a mock in the usual sense: it inherits the real client so the
    calling code under test is the calling code that ships, and only the
    transport is replaced.
    """

    def __init__(self, answers: Sequence[str], model: str = "scripted",
                 host: str = "scripted://none") -> None:
        LLMClient.__init__(self, model=model, host=host, backend=BACKEND_OLLAMA)
        self.answers = list(answers)
        #: Every prompt this client was asked, so a test can assert on what the
        #: model was allowed to see.
        self.prompts: List[Sequence[Dict[str, str]]] = []

    def detect_backend(self) -> Optional[str]:
        return BACKEND_OLLAMA

    def list_models(self) -> List[str]:
        return [self.model]

    def chat(self, messages: Sequence[Dict[str, str]], json_mode: bool = True,
             temperature: float = DEFAULT_TEMPERATURE) -> LLMResult:
        self.prompts.append(list(messages))
        if not self.answers:
            return LLMResult(ok=False, host=self.host, backend=BACKEND_OLLAMA,
                             model=self.model,
                             error="scripted client ran out of answers")
        text = self.answers.pop(0)
        return LLMResult(text=text, ok=bool(text), host=self.host,
                         backend=BACKEND_OLLAMA, model=self.model,
                         error="" if text else "the model returned an empty answer",
                         attempts=[text[:300]])


def client_from_config(config: Any) -> Optional[LLMClient]:
    """A client for this run, or None when the operator turned the model off.

    None is a legitimate configuration and not an error: ``--no-llm`` produces a
    rule-only report, which the design requires to be a supported mode rather
    than a degraded one.
    """
    if config is not None and not getattr(config, "use_llm", True):
        return None
    host = getattr(config, "llm_host", "") if config is not None else ""
    model = getattr(config, "llm_model", "") if config is not None else ""
    return LLMClient(model=model, host=host)
