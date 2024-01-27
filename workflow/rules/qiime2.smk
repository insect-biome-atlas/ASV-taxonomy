localrules:
    sintax2qiime_input,
    qiime2_import_seqs,

rule sintax2qiime_input:
    """
    Output should look like this:
    seqid1	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
    seqid2	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__

    fasta output should be:
    >seqid1
    ATGCGGGCTAGAGTAGCGAT...
    """
    input:
        tsv=lambda wildcards: config["qiime2"]["ref"][wildcards.ref]["taxonomy"],
        fasta=lambda wildcards: config["qiime2"]["ref"][wildcards.ref]["fasta"],
    output:
        tsv="results/qiime2/{ref}/taxonomy.tsv",
        fasta="results/qiime2/{ref}/seqs.fasta",
    log:
        "results/qiime2/{ref}/sintax2qiime_input.log"
    params:
        ranks = lambda wildcards: config["qiime2"]["ref"][wildcards.ref]["ranks"],
        output_dir = lambda wildcards, output: os.path.dirname(output.tsv),
    shell:
        """
        python workflow/scripts/sintax2qiime_input.py {input.tsv} {input.fasta} {output.tsv} {output.fasta} --ranks {params.ranks}> {log} 2>&1
        """
        

def qiime2_seqs(wildcards):
    if config["qiime2"]["ref"][wildcards.ref]["format"] == "sintax":
        return f"results/qiime2/{wildcards.ref}/seqs.fasta"
    else:
        return config["qiime2"]["ref"][wildcards.ref]["fasta"]

rule qiime2_import_seqs:
    output:
        "results/qiime2/{ref}/seqs.qza"
    input:
        qiime2_seqs,
    log:
        "results/qiime2/{ref}/qiime2_import_seqs.log"
    container:
        "docker://qiime2/core:2020.8"
    threads: 1
    resources:
        mem_mb = mem_allowed,
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Sequence]' --input-path {input} --output-path {output} > {log} 2>&1
        """

def qiime2_taxonomy(wildcards):
    if config["qiime2"]["ref"][wildcards.ref]["format"] == "sintax":
        return f"results/qiime2/{wildcards.ref}/taxonomy.tsv"
    else:
        return config["qiime2"]["ref"][wildcards.ref]["taxonomy"]

rule qiime2_import_taxonomy:
    output:
        "results/qiime2/{ref}/taxonomy.qza"
    input:
        qiime2_taxonomy,
    log:
        "results/qiime2/{ref}/qiime2_import_taxonomy.log"
    conda:
        "../envs/qiime2.yml"
    threads: 1
    resources:
        mem_mb = mem_allowed,
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat \
            --input-path {input} --output-path {output} > {log} 2>&1
        """

rule qiime2_train:
    output:
        "results/qiime2/{ref}/classifier.qza"
    input:
        tax="results/qiime2/{ref}/taxonomy.qza",
        seq="results/qiime2/{ref}/seqs.qza"
    log:
        "results/qiime2/{ref}/qiime2_train.log"
    conda:
        "../envs/qiime2.yml"
    resources:
        runtime=60 * 10,
        mem_mb = mem_allowed,
    threads: 1
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads {input.seq} --i-reference-taxonomy {input.tax} \
            --o-classifier {output} > {log} 2>&1
        """

rule qiime2_classify:
    output:
        "results/qiime2/{ref}/queries/{query}/taxonomy.qza"
    input:
        classifier="results/qiime2/{ref}/classifier.qza",
        qry=lambda wildcards: config["qiime2"]["query"][wildcards.query],
    log:
        "results/qiime2/{ref}/queries/{query}/qiime2_classify.log"
    threads: 1
    conda:
        "../envs/qiime2.yml"
    resources:
        runtime = 60*10,
        mem_mb = mem_allowed,
    shell:
        """
        qiime feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.qry} --o-classification {output} > {log} 2>&1
        """
