"""Alignment, variant calling, direct pileup and coverage — the measuring layer.

Nothing in here interprets anything. These modules turn reads into numbers:
a sorted BAM, a variant list, per-allele depths at named positions, and a
coverage summary. Every judgement about what those numbers mean belongs to
``resistance/``, ``typing/`` and ``contamination/``, which is why an engine
returns ``None`` for a quantity it could not measure rather than a zero that
reads like a measurement.
"""
