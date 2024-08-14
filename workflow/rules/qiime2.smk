localrules:
    sintax2qiime_input,
    qiime2_import_qry_seqs,
    qiime2_import_ref_seqs,
    qiime2_import_taxonomy,
    qiime2_export,
    split_qiime_input

wildcard_constraints:
    classifier = "vsearch|sklearn",

rule run_qiime2_vsearch:
    input:
        expand("results/qiime2/{ref}/queries/{query}/taxonomy_vsearch.tsv",
            ref=config["qiime2"]["ref"].keys(), query=config["qiime2"]["query"].keys())

rule run_qiime2_sklearn:
    input:
        expand("results/qiime2/{ref}/queries/{query}/taxonomy_sklearn.tsv",
            ref=config["qiime2"]["ref"].keys(), query=config["qiime2"]["query"].keys())

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

def qiime2_qry_seqs(wildcards):
    if config["qiime2"]["query"][wildcards.query]["format"] == "sintax":
        return f"results/qiime2/{wildcards.ref}/seqs.fasta"
    else:
        return config["qiime2"]["ref"][wildcards.ref]["fasta"]

def qiime2_ref_seqs(wildcards):
    if config["qiime2"]["ref"][wildcards.ref]["format"] == "sintax":
        return f"results/qiime2/{wildcards.ref}/seqs.fasta"
    else:
        return config["qiime2"]["ref"][wildcards.ref]["fasta"]


splits=[f'split{x:03d}' for x in list(range(1,1001))]

checkpoint split_qiime_input:
    """
    Splits the QIIME fasta file into fixed size chunks
    """
    output:
        directory("results/qiime2/{ref}/queries/{query}/splits")
    input:
        qry=lambda wildcards: config["qiime2"]["query"][wildcards.query],
    log:
        "logs/qiime2/qiime2.{ref}.{query}.split.log"
    params:
        outdir=lambda wildcards, output: output[0],
        size=500
    resources:
        runtime = 60,
    threads: 1
    shell:
        """
        cat {input.qry} | seqkit split2 -O {params.outdir} -j {threads} -s {params.size} >{log} 2>&1
        """

rule qiime2_import_ref_seqs:
    output:
        "results/qiime2/{ref}/seqs.qza"
    input:
        qiime2_ref_seqs,
    log:
        "results/qiime2/{ref}/qiime2_import_seqs.log"
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    resources:
        mem_mb = mem_allowed,
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Sequence]' --input-path {input} --output-path {output} > {log} 2>&1
        """
        
rule qiime2_import_qry_seqs:
    output:
        "results/qiime2/{ref}/queries/{query}/splits/{split}/{split}.qza"
    input:
        "results/qiime2/{ref}/queries/{query}/splits/stdin.part_{split}.fasta"
    log:
        "logs/qiime2/{ref}/queries/{query}/qiime2_import_seqs.{split}.log"
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
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
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    resources:
        mem_mb = mem_allowed,
        runtime = 60
    shell:
        """
        qiime tools import --type 'FeatureData[Taxonomy]' --input-format TSVTaxonomyFormat \
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
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    resources:
        runtime=60 * 24 * 10,
        mem_mb = mem_allowed,
    threads: 20
    shell:
        """
        qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads {input.seq} --i-reference-taxonomy {input.tax} \
            --o-classifier {output} > {log} 2>&1
        """

rule qiime2_classify_sklearn:
    output:
        "results/qiime2/{ref}/queries/{query}/splits/{split}/taxonomy_sklearn.qza"
    input:
        classifier="results/qiime2/{ref}/classifier.qza",
        qry=rules.qiime2_import_qry_seqs.output[0]
    log:
        "logs/qiime2/{ref}/queries/{query}/qiime2_classify_sklearn.{split}.log"
    threads: 20
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    resources:
        runtime = 60 * 10,
        mem_mb = mem_allowed,
    shell:
        """
        qiime feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.qry} --o-classification {output} > {log} 2>&1
        """

rule qiime2_classify_vsearch:
    output:
        vsearch="results/qiime2/{ref}/queries/{query}/splits/{split}/taxonomy_vsearch.qza",
        hits="results/qiime2/{ref}/queries/{query}/splits/{split}/taxonomy_hits.qza",
    input:
        ref="results/qiime2/{ref}/seqs.qza",
        ref_tax="results/qiime2/{ref}/taxonomy.qza",
        qry=rules.qiime2_import_qry_seqs.output[0]
    log:
        "logs/qiime2/{ref}/queries/{query}/qiime2_classify_vsearch.{split}.log"
    threads: 20
    container:
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    resources:
        runtime = 60 * 10,
        mem_mb = mem_allowed,
    shell:
        """
        qiime feature-classifier classify-consensus-vsearch --i-reference-reads {input.ref} --i-query {input.qry} \
            --i-reference-taxonomy {input.ref_tax} --o-classification {output.vsearch} --o-search-results {output.hits} \
            --p-threads {threads} --verbose > {log} 2>&1
        """

rule qiime2_export:
    output:
        "results/qiime2/{ref}/queries/{query}/splits/{split}/taxonomy_{classifier}.tsv"
    input:
        "results/qiime2/{ref}/queries/{query}/splits/{split}/taxonomy_{classifier}.qza"
    log:
        "logs/qiime2/{ref}/queries/{query}/qiime2_export_{classifier}.{split}.log"
    container: 
        "docker://quay.io/qiime2/core:2023.9" # may have to be built as part of a SLURM job on Uppmax
    threads: 1
    shell:
        """
        qiime tools export --input-path {input} --output-path {output[0]} --output-format TSVTaxonomyFormat > {log} 2>&1
        """

def aggregate_qiime(wildcards):
    checkpoint_output = checkpoints.split_qiime_input.get(**wildcards).output[0]
    return expand("results/qiime2/{ref}/queries/{query}/splits/{split}/taxonomy_{classifier}.tsv",
                    ref=wildcards.ref, query=wildcards.query, classifier=wildcards.classifier, 
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)

rule collate_qiime:
    """
    Concatenates the qiime output files into a single file
    """
    output:
        "results/qiime2/{ref}/queries/{query}/taxonomy_{classifier}.tsv"
    input:
        aggregate_qiime,
    run:
        with open(output[0], "w") as out:
            for i, f in enumerate(input):
                with open(f, 'r') as infile:
                    for j, line in enumerate(infile):
                        if j == 0 and i == 0:
                            out.write(line)
                        elif j > 0:
                            out.write(line)

rule parse_qiime:
    output:
        "results/qiime2/{ref}/queries/{query}/taxonomy_{classifier}_parsed.tsv"
    input:
        rules.collate_qiime.output[0]
    run:
        import pandas as pd
        from workflow.scripts.evaluate_classifier import parse_qiime2
        df = pd.read_csv(input[0], sep="\t", index_col=0)
        parsed = parse_qiime2(df)
        parsed.to_csv(output[0], sep="\t")
    