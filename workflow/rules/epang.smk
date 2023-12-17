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
    params:
        src=srcdir("../scripts/nexus2newick.py")
    envmodules:
        "bioinfo-tools",
        "biopython"
    shell:
        """
        python {params.src} {input} {output} >{log} 2>&1
        """

def ref_tree(wildcards):
    if config["epa-ng"]["ref"][wildcards.ref]["tree_format"] == "nexus":
        return rules.nexus2newick.output[0]
    elif config["epa-ng"]["ref"][wildcards.ref]["tree_format"] == "newick":
        return config["epa-ng"]["ref"][wildcards.ref]["tree"]


rule extract_ref_taxonomy:
    output:
        "results/epa-ng/{ref}/taxon_file.tsv",
    input:
        ref_tree,
    log:
        "logs/epa-ng/{ref}/extract_ref_taxonomy.log",
    params:
        ranks=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["tree_ranks"],
        src=srcdir("../scripts/extract_ref_taxonomy.py")
    shell:
        """
        python {params.src} {input} {output} --ranks {params.ranks} >{log} 2>&1
        """


rule nexus2fasta:
    output:
        "results/epa-ng/{ref}/{ref}.fasta",
    input:
        lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["msa"],
    log:
        "logs/epa-ng/{ref}/nexus2fasta.log",
    params:
        src=srcdir("../scripts/convertalign.py")
    envmodules:
        "bioinfo-tools",
        "biopython"
    shell:
        """
        python {params.src} {input} nexus {output} fasta >{log} 2>&1
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
        ref_msa=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["msa"],
    log:
        "logs/epa-ng/{ref}/hmmalign.{query}.log",
    resources:
        runtime=120,
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
        ref_msa="results/epa-ng/{ref}/hmmalign/{query}/ref.fasta",
        qry_msa="results/epa-ng/{ref}/hmmalign/{query}/query.fasta",
    input:
        ref_msa=ref_msa,
        msa=rules.hmm_align.output,
    log:
        "logs/epa-ng/{ref}/split_aln.{query}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.ref_msa),
    shell:
        """
        epa-ng --redo --out-dir {params.outdir} --split {input.ref_msa} {input.msa} > {log} 2>&1
        mv {params.outdir}/query.fasta {output.qry_msa}
        mv {params.outdir}/reference.fasta {output.ref_msa}
        """

rule raxml_evaluate:
    output:
        "results/epa-ng/{ref}/raxml-ng/{query}/info.raxml.bestModel",
    input:
        tree=rules.nexus2newick.output,
        msa=rules.split_aln.output.ref_msa,
    log:
        "logs/epa-ng/{ref}/raxml-ng.{query}.log",
    params:
        model=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["model"],
        prefix=lambda wildcards, output: os.path.dirname(output[0]) + "/info",
    threads: 4
    resources:
        runtime=60,
    shell:
        """
        raxml-ng --threads {threads} --evaluate --msa {input.msa} --tree {input.tree} --prefix {params.prefix} --model {params.model} >{log} 2>&1
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
    threads: 4
    resources:
        runtime=60,
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
    if config["epa-ng"]["gappa"]["distribution-ratio"] == -1:
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
        consensus_thresh=config["epa-ng"]["gappa"]["consensus-thresh"],
        distribution_ratio=get_dist_ratio(config),
    threads: 4
    resources:
        runtime=120,
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
        src=srcdir("../scripts/gappa2taxdf.py"),
    shell:
        """
        python {params.src} {input} {output} --ranks {params.ranks} >{log} 2>&1
        """


rule write_config:
    output:
        "results/epa-ng/{ref}/queries/{query}/config.yml",
    run:
        import yaml
        with open(output[0], "w") as fhout:
            yaml.safe_dump(config, fhout, default_flow_style=False, sort_keys=False)


rule write_software:
    output:
        "results/epa-ng/{ref}/queries/{query}/software.txt",
    shell:
        """
        conda list > {output}
        """