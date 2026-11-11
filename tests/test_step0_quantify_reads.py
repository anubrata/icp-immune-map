"""Raw reads -> gene counts: adding Salmon's transcript estimates up per gene, and reading GENCODE files."""

import gzip

import numpy as np
import pandas as pd

from step0_build_salmon_index import gene_table, sequence_names, transcript_table
from step0_quantify_reads import gene_counts

TX2GENE = pd.DataFrame({"transcript_id": ["A-1", "A-2", "B-1", "C-1"], "gene_id": ["A", "A", "B", "C"]})


def quant(tpm, length, reads):
    """A quant.sf-like table for transcripts A-1, A-2, B-1, C-1."""
    return pd.DataFrame({"TPM": tpm, "EffectiveLength": length, "NumReads": reads}, index=["A-1", "A-2", "B-1", "C-1"])


def test_counts_match_reads_when_lengths_are_stable():
    """With the same transcript lengths everywhere, lengthScaledTPM gives back each gene's summed reads."""
    length = [1000.0, 500.0, 2000.0, 800.0]
    reads = {"s1": [300.0, 100.0, 800.0, 0.0], "s2": [150.0, 50.0, 1600.0, 0.0]}
    quants = {}
    for sample, r in reads.items():
        rate = np.array(r) / np.array(length)               # reads per base ~ abundance
        quants[sample] = quant(rate / rate.sum() * 1e6, length, r)
    counts = gene_counts(quants, TX2GENE)
    assert counts.dtypes.eq("int64").all()
    assert counts.loc["A"].tolist() == [400, 200] and counts.loc["B"].tolist() == [800, 1600]
    assert counts.loc["C"].tolist() == [0, 0]                # never expressed: 0, not missing
    assert (counts.sum() == pd.Series({"s1": 1200, "s2": 1800})).all()


def test_isoform_switch_does_not_look_like_a_change():
    """Gene A has the same abundance (TPM) in both samples, but s2 uses its short isoform. Summed reads would
    call A lower in s2; lengthScaledTPM uses A's average length, so A gets the same share of counts."""
    q1 = quant([100.0, 0.0, 100.0, 0.0], [1000.0, 500.0, 1000.0, 800.0], [1000.0, 0.0, 1000.0, 0.0])
    q2 = quant([0.0, 100.0, 100.0, 0.0], [1000.0, 500.0, 1000.0, 800.0], [0.0, 500.0, 1000.0, 0.0])
    counts = gene_counts({"s1": q1, "s2": q2}, TX2GENE)
    share = counts.loc["A"] / counts.sum()
    assert abs(share["s1"] - share["s2"]) < 1e-3


def test_reads_gencode_files(tmp_path):
    fasta = tmp_path / "transcripts.fa.gz"
    with gzip.open(fasta, "wt") as f:
        f.write(">ENST01.2|ENSG01.5|-|-|GENE1-201|GENE1|1657|lncRNA|\nACGT\n>ENST02.1|ENSG02.1|-|-|GENE2-201|GENE2|900|protein_coding|\nAC\n")
    assert transcript_table(fasta).values.tolist() == [["ENST01.2", "ENSG01.5"], ["ENST02.1", "ENSG02.1"]]
    assert sequence_names(fasta) == ["ENST01.2|ENSG01.5|-|-|GENE1-201|GENE1|1657|lncRNA|", "ENST02.1|ENSG02.1|-|-|GENE2-201|GENE2|900|protein_coding|"]

    gtf = tmp_path / "annotation.gtf.gz"
    with gzip.open(gtf, "wt") as f:
        f.write("##description: test\n"
                'chr17\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG01.5"; gene_type "protein_coding"; gene_name "TP53";\n'
                'chr17\tHAVANA\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG01.5"; transcript_id "ENST01.2";\n'
                'chrM\tENSEMBL\tgene\t1\t50\t.\t+\t.\tgene_id "ENSG02.1"; gene_type "Mt_tRNA"; gene_name "MT-TF";\n')
    assert gene_table(gtf).values.tolist() == [["ENSG01.5", "TP53", "protein_coding", "17"], ["ENSG02.1", "MT-TF", "Mt_tRNA", "MT"]]


def test_source_setting(tmp_path):
    """Datasets not on GEO name their source; it must have a label and a link."""
    import pytest
    from common import ConfigError, load_config
    from conftest import write_config
    source = {"label": "BioProject PRJNA846869", "url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA846869"}
    assert load_config(write_config(tmp_path / "a.yaml", tmp_path / "x.txt", source=source))["source"] == source
    assert load_config(write_config(tmp_path / "b.yaml", tmp_path / "x.txt"))["source"] is None
    with pytest.raises(ConfigError, match="source needs"):
        load_config(write_config(tmp_path / "c.yaml", tmp_path / "x.txt", source={"label": "no link"}))
