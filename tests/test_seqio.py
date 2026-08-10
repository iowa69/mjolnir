"""Input detection: pairing, sample naming and platform, on real filenames.

The filenames in these tests are the ones on the drives this tool was written
for. ``30-20_S1_R1_001.fastq.gz`` and ``226-18_S8_L001_R1_001.fastq.gz`` are
bcl2fastq output, and both are named wrongly by the obvious rules: "strip after
the last underscore" gives ``30-20_S1_R1``, and "split on the first underscore"
gives ``30``. A sample renamed at the input stage is renamed in the report, in
the joint variant table and in the cluster, so this is not cosmetic.

The second rule tested here is that a read number is only believed when its
partner is present. A lone ``patient_2.fastq.gz`` is a single-end file for a
sample called ``patient_2``, not read 2 of a pair — the other reading silently
renames the sample to ``patient``.
"""

from __future__ import annotations

import gzip

import pytest

from mjolnir import seqio
from mjolnir.utils import MjolnirError

#: One tiny well-formed FASTQ record. Long enough to look like Illumina, short
#: enough that nothing here touches a real read file.
READ = ("@M03605:1:000000000-ABCDE:1:1101:1:1 1:N:0:1\n"
        + "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        + "+\n" + ("I" * 100) + "\n")


def _fastq(tmp_path, name, records=4):
    path = tmp_path / name
    payload = (READ * records).encode()
    if name.endswith(".gz"):
        with gzip.open(str(path), "wb") as handle:
            handle.write(payload)
    else:
        path.write_bytes(payload)
    return path


# ------------------------------------------------------------- sample naming

@pytest.mark.parametrize("name,expected", [
    ("30-20_S1_R1_001.fastq.gz", "30-20"),
    ("30-20_S1_R2_001.fastq.gz", "30-20"),
    ("226-18_S8_L001_R1_001.fastq.gz", "226-18"),
    ("226-18_S8_L001_R2_001.fastq.gz", "226-18"),
    ("M_chimaera_TN_S12_L001_R1_001.fastq.gz", "M_chimaera_TN"),
    ("sample_R1.fastq.gz", "sample"),
    ("sample.R2.fq.gz", "sample"),
])
def test_sample_name_from_real_filenames(name, expected):
    assert seqio.sample_name(name) == expected


def test_the_sample_sheet_index_is_decoration_not_part_of_the_name():
    parsed = seqio.parse_read_name("226-18_S8_L001_R1_001.fastq.gz")
    assert parsed.sample == "226-18"
    assert parsed.read == 1
    assert parsed.lane == "001"
    assert parsed.convention == seqio.CONVENTION_BCL2FASTQ


def test_a_bare_numeric_suffix_is_not_stripped_on_its_own_evidence():
    """``patient_2`` is a plausible isolate name; renaming it would be silent."""
    assert seqio.sample_name("patient_2.fastq.gz") == "patient_2"
    assert seqio.parse_read_name("patient_2.fastq.gz").ambiguous is True


# ----------------------------------------------------------------- pairing

def test_the_local_illumina_filenames_pair(tmp_path):
    paths = [_fastq(tmp_path, "30-20_S1_R1_001.fastq.gz"),
             _fastq(tmp_path, "30-20_S1_R2_001.fastq.gz")]
    groups = seqio.group_reads(paths)
    assert len(groups) == 1
    group = groups[0]
    assert group.sample_id == "30-20"
    assert group.paired is True
    assert [p.name for p in group.files] == ["30-20_S1_R1_001.fastq.gz",
                                             "30-20_S1_R2_001.fastq.gz"]
    assert group.convention == seqio.CONVENTION_BCL2FASTQ


def test_the_lane_split_filenames_pair_too(tmp_path):
    paths = [_fastq(tmp_path, "226-18_S8_L001_R1_001.fastq.gz"),
             _fastq(tmp_path, "226-18_S8_L001_R2_001.fastq.gz")]
    groups = seqio.group_reads(paths)
    assert [g.sample_id for g in groups] == ["226-18"]
    assert groups[0].paired is True


def test_two_samples_in_one_directory_stay_two_samples(tmp_path):
    paths = [_fastq(tmp_path, "30-20_S1_R1_001.fastq.gz"),
             _fastq(tmp_path, "30-20_S1_R2_001.fastq.gz"),
             _fastq(tmp_path, "226-18_S8_L001_R1_001.fastq.gz"),
             _fastq(tmp_path, "226-18_S8_L001_R2_001.fastq.gz")]
    groups = seqio.group_reads(paths)
    assert sorted(g.sample_id for g in groups) == ["226-18", "30-20"]
    assert all(g.paired for g in groups)


def test_a_sample_split_across_lanes_is_refused_with_the_command_to_fix_it(tmp_path):
    paths = [_fastq(tmp_path, "226-18_S8_L001_R1_001.fastq.gz"),
             _fastq(tmp_path, "226-18_S8_L002_R1_001.fastq.gz")]
    with pytest.raises(MjolnirError) as excinfo:
        seqio.group_reads(paths)
    message = str(excinfo.value)
    assert "split across lanes" in message
    assert "cat " in message, "the refusal has to say how to proceed"
    assert "duplicate marking" in message


def test_a_lone_r1_is_single_end_under_its_own_name(tmp_path):
    groups = seqio.group_reads([_fastq(tmp_path, "30-20_S1_R1_001.fastq.gz")])
    assert len(groups) == 1
    assert groups[0].sample_id == "30-20"
    assert groups[0].paired is False
    assert any("mate was not given" in note for note in groups[0].notes)


def test_a_lone_bare_suffix_file_keeps_its_whole_name(tmp_path):
    groups = seqio.group_reads([_fastq(tmp_path, "patient_2.fastq.gz")])
    assert groups[0].sample_id == "patient_2"
    assert any("rather than as read 2" in note for note in groups[0].notes)


def test_a_paired_bare_suffix_set_is_believed(tmp_path):
    paths = [_fastq(tmp_path, "patient_1.fastq.gz"),
             _fastq(tmp_path, "patient_2.fastq.gz")]
    groups = seqio.group_reads(paths)
    assert [g.sample_id for g in groups] == ["patient"]
    assert groups[0].paired is True


def test_two_files_for_the_same_read_are_refused(tmp_path):
    paths = [_fastq(tmp_path, "s_R1.fastq.gz"), _fastq(tmp_path, "s.R1.fastq.gz")]
    with pytest.raises(MjolnirError) as excinfo:
        seqio.group_reads(paths)
    assert "more than one file for read 1" in str(excinfo.value)


def test_index_reads_are_dropped_and_the_drop_is_stated(tmp_path):
    paths = [_fastq(tmp_path, "30-20_S1_R1_001.fastq.gz"),
             _fastq(tmp_path, "30-20_S1_R2_001.fastq.gz"),
             _fastq(tmp_path, "30-20_S1_I1_001.fastq.gz")]
    groups = seqio.group_reads(paths)
    assert len(groups) == 1
    assert len(groups[0].files) == 2
    assert any("index read" in note for note in groups[0].notes)


# ------------------------------------------------------------ format sniffing

def test_a_gzipped_fastq_is_recognised(tmp_path):
    path = _fastq(tmp_path, "30-20_S1_R1_001.fastq.gz")
    assert seqio.sniff_format(path) == seqio.FORMAT_FASTQ
    assert seqio.is_fastq(path) is True
    assert seqio.compression(path) == "gzip"


def test_a_fasta_is_recognised(tmp_path):
    path = tmp_path / "assembly.fasta"
    path.write_text(">contig_1\nACGTACGTACGT\n")
    assert seqio.sniff_format(path) == seqio.FORMAT_FASTA
    assert seqio.is_fasta(path) is True


def test_a_file_whose_extension_lies_about_compression_is_caught(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(READ.encode())  # named .gz, not actually gzipped
    assert seqio.compression_mismatch(path)


def test_an_empty_fastq_is_not_a_valid_input(tmp_path):
    path = tmp_path / "empty.fastq"
    path.write_text("")
    with pytest.raises(MjolnirError):
        seqio.validate_fastq(path)


# ----------------------------------------------------------- platform detection

def test_short_reads_are_called_illumina_from_the_reads_not_the_name(tmp_path):
    paths = [_fastq(tmp_path, "30-20_S1_R1_001.fastq.gz", records=200),
             _fastq(tmp_path, "30-20_S1_R2_001.fastq.gz", records=200)]
    evidence = seqio.detect_platform(paths, paired=True)
    assert evidence.platform == "illumina"
    assert evidence.reasons, "the call must say what it was based on"


def test_an_assembly_is_detected_as_fasta(tmp_path):
    path = tmp_path / "M_chimaera_TN.fasta"
    path.write_text(">contig_1\n" + ("ACGT" * 500) + "\n")
    evidence = seqio.detect_platform([path])
    assert evidence.platform == "fasta"
