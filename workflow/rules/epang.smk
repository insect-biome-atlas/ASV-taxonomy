localrules:
    nexus2newick,
    extract_ref_taxonomy,
    nexus2fasta,
    split_aln,
    gappa2taxdf,
    write_config,
    write_software,

rule run_epa_ng:
    input:
        expand("results/epa-ng/{ref}/queries/{query}/{f}",
            ref = config["epa-ng"]["ref"].keys(), query = config["epa-ng"]["query"].keys(),
            f = ["software.txt", "config.yml"]),
        expand("results/epa-ng/{ref}/queries/{query}.taxonomy.tsv",
            ref = config["epa-ng"]["ref"].keys(), query = config["epa-ng"]["query"].keys())

rule nexus2newick:
    output:
        "results/epa-ng/{ref}/tree.nwk",
    input:
        lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["tree"],
    log:
        "logs/epa-ng/{ref}/nexus2newick.log",
    envmodules:
        "bioinfo-tools",
        "biopython"
    shell:
        """
        python workflow/scripts/nexus2newick.py {input} {output} >{log} 2>&1
        """

def ref_tree(wildcards):
    if config["epa-ng"]["ref"][wildcards.ref]["tree_format"] == "nexus":
        return rules.nexus2newick.output[0]
    elif config["epa-ng"]["ref"][wildcards.ref]["tree_format"] == "newick":
        return config["epa-ng"]["ref"][wildcards.ref]["tree"]

def ref_msa(wildcards):
    if config["epa-ng"]["ref"][wildcards.ref]["msa_format"] == "nexus":
        return rules.nexus2fasta.output[0]
    elif config["epa-ng"]["ref"]["wildcards.ref"]["msa_format"] == "fasta":
        return config["epa-ng"]["ref"][wildcards.ref]["msa"]

rule extract_ref_taxonomy:
    output:
        "results/epa-ng/{ref}/taxon_file.tsv",
    input:
        ref_tree,
    log:
        "logs/epa-ng/{ref}/extract_ref_taxonomy.log",
    params:
        ranks=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["tree_ranks"],
    shell:
        """
        python workflow/scripts/extract_ref_taxonomy.py {input} {output} --ranks {params.ranks} >{log} 2>&1
        """


rule nexus2fasta:
    output:
        "results/epa-ng/{ref}/{ref}.fasta",
    input:
        lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["msa"],
    log:
        "logs/epa-ng/{ref}/nexus2fasta.log",
    envmodules:
        "bioinfo-tools",
        "biopython"
    shell:
        """
        python workflow/scripts/convertalign.py {input} nexus {output} fasta >{log} 2>&1
        """

def ref_msa(wildcards):
    if config["epa-ng"]["ref"][wildcards.ref]["msa_format"] == "nexus":
        return rules.nexus2fasta.output[0]
    elif config["epa-ng"]["ref"][wildcards.ref]["msa_format"] == "fasta":
        return config["epa-ng"]["ref"][wildcards.ref]["msa"]

rule hmm_build:
    output:
        "results/epa-ng/{ref}/hmmalign/{ref}.hmm",
    input:
        ref_msa,
    log:
        "logs/epa-ng/{ref}/hmmbuild.log",
    resources:
        runtime=60,
    envmodules:
        "bioinfo-tools",
        "hmmer/3.3.2"
    shell:
        """
        hmmbuild {output} {input} > {log} 2>&1
        """

rule hmm_align:
    output:
        "results/epa-ng/{ref}/hmmalign/{query}/{ref}.{query}.fasta",
    input:
        hmm=rules.hmm_build.output,
        qry=lambda wildcards: config["epa-ng"]["query"][wildcards.query],
        ref_msa=ref_msa
    log:
        "logs/epa-ng/{ref}/hmmalign.{query}.log",
    resources:
        runtime=60*24*10,
        mem_mb=mem_allowed,
    envmodules:
        "bioinfo-tools",
        "hmmer/3.3.2"
    threads: 4
    shell:
        """
        hmmalign --trim --mapali {input.ref_msa} --outformat afa -o {output} {input.hmm} {input.qry} > {log} 2>&1
        """

rule split_aln:
    output:
        ref_msa="results/epa-ng/{ref}/hmmalign/{query}/reference.fasta",
        qry_msa="results/epa-ng/{ref}/hmmalign/{query}/query.fasta",
    input:
        ref_msa=ref_msa,
        msa=rules.hmm_align.output,
    log:
        "logs/epa-ng/{ref}/split_aln.{query}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.ref_msa),
    envmodules:
        "bioinfo-tools",
        "EPA-ng/0.3.8"
    shell:
        """
        epa-ng --redo --out-dir {params.outdir} --split {input.ref_msa} {input.msa} > {log} 2>&1
        """

rule raxml_evaluate:
    output:
        "results/epa-ng/{ref}/raxml-ng/{query}/info.raxml.bestModel",
    input:
        tree=ref_tree,
        msa=rules.split_aln.output.ref_msa,
    log:
        "logs/epa-ng/{ref}/raxml-ng.{query}.log",
    params:
        model=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["model"],
        prefix=lambda wildcards, output: os.path.dirname(output[0]) + "/info",
    envmodules:
        "bioinfo-tools",
        "RAxML-NG/1.1.0"
    threads: 2
    resources:
        runtime=60,
        mem_mb=mem_allowed,
    shell:
        """
        raxml-ng --redo --threads {threads} --evaluate --msa {input.msa} --tree {input.tree} --prefix {params.prefix} --model {params.model} >{log} 2>&1
        """

rule epa_ng:
    output:
        "results/epa-ng/{ref}/queries/{query}/epa_result.jplace",
    input:
        qry=rules.split_aln.output.qry_msa,
        ref_msa=rules.split_aln.output.ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output,
    log:
        "logs/epa-ng/{ref}/epa-ng.{query}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    envmodules:
        "bioinfo-tools",
        "EPA-ng/0.3.8"
    threads: 20
    resources:
        runtime=60*24,
        mem_mb=mem_allowed,
    shell:
        """
        epa-ng --redo -T {threads} --tree {input.ref_tree} --ref-msa {input.ref_msa} \
            --query {input.qry} --out-dir {params.outdir} --model {input.info} >{log} 2>&1
        """


def ref_taxonomy(wildcards):
    if config["epa-ng"]["ref"][wildcards.ref]["ref_taxonomy"]:
        return config["epa-ng"]["ref"][wildcards.ref]["ref_taxonomy"]
    else:
        return rules.extract_ref_taxonomy.output


def get_dist_ratio(config):
    if config["epa-ng"]["gappa"]["distribution_ratio"] == -1:
        return ""
    else:
        return f"--distribution-ratio {config['gappa']['distribution_ratio']}"

rule gappa_assign:
    output:
        "results/epa-ng/{ref}/queries/{query}/per_query.tsv",
    input:
        json=rules.epa_ng.output,
        taxonfile=ref_taxonomy,
    log:
        "logs/epa-ng/{ref}/gappa_assign.{query}.log",
    params:
        ranks_string=lambda wildcards: "|".join(config["epa-ng"]["ref"][wildcards.ref]["tree_ranks"]),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        consensus_thresh=config["epa-ng"]["gappa"]["consensus_thresh"],
        distribution_ratio=get_dist_ratio(config),
    #envmodules:
    #    "bioinfo-tools",
    #    "gappa/0.7.1"
    conda:
        "../envs/gappa.yml"
    threads: 20
    resources:
        runtime=60 *2,
        mem_mb=mem_allowed,
    shell:
        """
        gappa examine assign --threads {threads} --out-dir {params.outdir} \
            --jplace-path {input.json} --taxon-file {input.taxonfile} \
            --ranks-string '{params.ranks_string}' --per-query-results \
            --consensus-thresh {params.consensus_thresh} {params.distribution_ratio} \
            --best-hit --allow-file-overwriting > {log} 2>&1
        """


rule gappa2taxdf:
    output:
        "results/epa-ng/{ref}/queries/{query}.taxonomy.tsv",
    input:
        rules.gappa_assign.output[0],
    log:
        "logs/epa-ng/{ref}/gappa2taxdf.{query}.log"
    params:
        ranks=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["tree_ranks"],
    shell:
        """
        python workflow/scripts/gappa2taxdf.py {input} {output} --ranks {params.ranks} >{log} 2>&1
        """
