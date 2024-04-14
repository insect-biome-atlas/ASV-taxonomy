localrules:
    nexus2newick,
    extract_ref_taxonomy,
    nexus2fasta,
    split_aln,
    gappa2taxdf,
    write_config,
    write_software,

wildcard_constraints:
    heur="baseball-heur|dyn-heur|no-heur",
    placer="pplacer|epa-ng",

## target rules

def get_phylo_input(wildcards):
    input = []
    for ref in config["phylogeny"]["ref"].keys():
        for placement_tool in config["phylogeny"]["placement-tools"]:
            for heur in config["phylogeny"][placement_tool]["heuristics"]:
                for query in config["phylogeny"]["query"].keys():
                    input.append(f"results/{placement_tool}/{ref}/queries/{query}/{heur}/taxonomy.tsv")
    return input

rule run_phylo:
    input:
        get_phylo_input,

rule run_epa_ng:
    input:
        expand("results/epa-ng/{ref}/queries/{query}/{heur}/taxonomy.tsv",
            ref = config["phylogeny"]["ref"].keys(), 
            query = config["phylogeny"]["query"].keys(), 
            heur = config["phylogeny"]["epa-ng"]["heuristics"]
            )

rule run_pplacer:
    input:
        expand("results/pplacer/{ref}/queries/{query}/{heur}/taxonomy.tsv",
            ref = config["phylogeny"]["ref"].keys(), 
            query = config["phylogeny"]["query"].keys(), 
            heur = config["phylogeny"]["pplacer"]["heuristics"],
            )

rule nexus2newick:
    output:
        "results/phylogeny/{ref}/tree.nwk",
    input:
        lambda wildcards: config["phylogeny"]["ref"][wildcards.ref]["tree"],
    log:
        "logs/phylogeny/{ref}/nexus2newick.log",
    envmodules:
        "bioinfo-tools",
        "biopython"
    shell:
        """
        python workflow/scripts/nexus2newick.py {input} {output} >{log} 2>&1
        """

def ref_tree(wildcards):
    if config["phylogeny"]["ref"][wildcards.ref]["tree_format"] == "nexus":
        return rules.nexus2newick.output[0]
    elif config["phylogeny"]["ref"][wildcards.ref]["tree_format"] == "newick":
        return config["phylogeny"]["ref"][wildcards.ref]["tree"]

def ref_msa(wildcards):
    if config["phylogeny"]["ref"][wildcards.ref]["msa_format"] == "nexus":
        return rules.nexus2fasta.output[0]
    elif config["phylogeny"]["ref"]["wildcards.ref"]["msa_format"] == "fasta":
        return config["phylogeny"]["ref"][wildcards.ref]["msa"]

rule extract_ref_taxonomy:
    output:
        "results/phylogeny/{ref}/taxon_file.tsv",
    input:
        ref_tree,
    log:
        "logs/phylogeny/{ref}/extract_ref_taxonomy.log",
    params:
        ranks=lambda wildcards: config["phylogeny"]["ref"][wildcards.ref]["tree_ranks"],
    shell:
        """
        python workflow/scripts/extract_ref_taxonomy.py {input} {output} --ranks {params.ranks} >{log} 2>&1
        """

rule nexus2fasta:
    output:
        "results/phylogeny/{ref}/{ref}.fasta",
    input:
        lambda wildcards: config["phylogeny"]["ref"][wildcards.ref]["msa"],
    log:
        "logs/phylogeny/{ref}/nexus2fasta.log",
    envmodules:
        "bioinfo-tools",
        "biopython"
    shell:
        """
        python workflow/scripts/convertalign.py {input} nexus {output} fasta >{log} 2>&1
        """

def ref_msa(wildcards):
    if config["phylogeny"]["ref"][wildcards.ref]["msa_format"] == "nexus":
        return rules.nexus2fasta.output[0]
    elif config["phylogeny"]["ref"][wildcards.ref]["msa_format"] == "fasta":
        return config["phylogeny"]["ref"][wildcards.ref]["msa"]

rule hmm_build:
    output:
        "results/phylogeny/{ref}/hmmalign/{ref}.hmm",
    input:
        ref_msa,
    log:
        "logs/phylogeny/{ref}/hmmbuild.log",
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
        "results/phylogeny/{ref}/hmmalign/{query}/{ref}.{query}.fasta",
    input:
        hmm=rules.hmm_build.output,
        qry=lambda wildcards: config["phylogeny"]["query"][wildcards.query],
        ref_msa=ref_msa
    log:
        "logs/phylogeny/{ref}/hmmalign.{query}.log",
    resources:
        runtime=60*24*10,
        mem_mb=mem_allowed,
    params:
        tmpdir=lambda wildcards: f"$TMPDIR/{wildcards.ref}.{wildcards.query}.raxml-ng",
    envmodules:
        "bioinfo-tools",
        "hmmer/3.3.2"
    threads: 4
    shell:
        """
        mkdir -p {params.tmpdir}
        hmmalign --trim --mapali {input.ref_msa} --outformat afa -o {params.tmpdir}/output.msa {input.hmm} {input.qry} > {log} 2>&1
        cut -f1 -d ' ' {params.tmpdir}/output.msa > {output[0]}
        rm -r {params.tmpdir}
        """

rule split_aln:
    output:
        ref_msa="results/phylogeny/{ref}/hmmalign/{query}/reference.fasta",
        qry_msa="results/phylogeny/{ref}/hmmalign/{query}/query.fasta",
    input:
        ref_msa=ref_msa,
        msa=rules.hmm_align.output,
    log:
        "logs/phylogeny/{ref}/split_aln.{query}.log",
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
        "results/phylogeny/{ref}/raxml-ng/{query}/info.raxml.bestModel",
    input:
        tree=ref_tree,
        msa=rules.split_aln.output.ref_msa,
    log:
        "logs/phylogeny/{ref}/raxml-ng.{query}.log",
    params:
        model=lambda wildcards: config["phylogeny"]["ref"][wildcards.ref]["model"],
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

## epa-ng

def get_heuristic(wildcards):
    if wildcards.heur == "no-heur":
        return "--no-heur"
    elif wildcards.heur == "dyn-heur":
        return "--dyn-heur 0.99999"
    elif wildcards.heur == "baseball-heur":
        return "--baseball-heur"

rule epa_ng:
    output:
        "results/epa-ng/{ref}/queries/{query}/{heur}/epa-ng_result.jplace",
    input:
        qry=rules.split_aln.output.qry_msa,
        ref_msa=rules.split_aln.output.ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output[0],
    log:
        "logs/epa-ng/{ref}/queries/{query}/{heur}/epa-ng.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        heur=get_heuristic,
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
            --query {input.qry} --out-dir {params.outdir} {params.heur} --model {input.info} >{log} 2>&1
        mv {params.outdir}/epa_result.jplace {output[0]}
        """

## pplacer

rule taxit_create:
    """
    Create refpkg for the tree with taxit, for later use with pplacer
    """
    output:
        "results/pplacer/{ref}/queries/{query}/taxit/refpkg/CONTENTS.json"
    input:
        ref_msa=ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output[0],
    log:
        "logs/taxit/{ref}/queries/{query}/create.log"
    conda: "../envs/pplacer.yml"
    #envmodules:
    #    "bioinfo-tools",
    #    "pplacer/1.1.alpha19"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
    shell:
        """
        taxit create -l {wildcards.ref} -P {params.outdir} --aln-fasta {input.ref_msa} --tree-stats {input.info} -tree-file {input.ref_tree} 2> {log}
        """

rule pplacer:
    """
    Run pplacer on query sequences against reference tree
    """
    output:
        jplace="results/pplacer/{ref}/queries/{query}/{heur}/pplacer_result.jplace"
    input:
        refpkg=rules.taxit_create.output[0],
        msa=rules.hmm_align.output[0]
    log:
        "logs/pplacer/{ref}/{query}/{heur}.log"
    envmodules:
        "bioinfo-tools",
        "pplacer/1.1.alpha19"
    params:
        max_strikes = lambda wildcards: 6 if wildcards.heur == "baseball" else 0
    resources:
        runtime=60*24*10,
        mem_mb=mem_allowed,
    shell:
        """
        pplacer -c {input.refpkg} {input.msa} -o {output.jplace} --timing --max-strikes {params.max_strikes} > {log} 2>&1
        """

## gappa

def ref_taxonomy(wildcards):
    if config["phylogeny"]["ref"][wildcards.ref]["ref_taxonomy"]:
        return config["phylogeny"]["ref"][wildcards.ref]["ref_taxonomy"]
    else:
        return rules.extract_ref_taxonomy.output


def get_dist_ratio(config):
    if config["phylogeny"]["gappa"]["distribution_ratio"] == -1:
        return ""
    else:
        return f"--distribution-ratio {config['phylogeny']['gappa']['distribution_ratio']}"

rule gappa_assign:
    """
    Run gappa taxonomic assignment on placement file
    Here {placer} is a wildcard that can be either 'pplacer' or 'epa-ng'
    """
    output:
        "results/{placer}/{ref}/queries/{query}/{heur}/per_query.tsv",
        "results/{placer}/{ref}/queries/{query}/{heur}/profile.tsv",
        "results/{placer}/{ref}/queries/{query}/{heur}/labelled_tree.newick",
    input:
        jplace="results/{placer}/{ref}/queries/{query}/{heur}/{placer}_result.jplace",
        taxonfile=ref_taxonomy,
    log:
        "logs/{placer}/{ref}/queries/{query}/{heur}/gappa_assign.log",
    params:
        ranks_string=lambda wildcards: "|".join(config["phylogeny"]["ref"][wildcards.ref]["tree_ranks"]),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        consensus_thresh=config["phylogeny"]["gappa"]["consensus_thresh"],
        distribution_ratio=get_dist_ratio(config),
    conda:
        "../envs/gappa.yml"
    threads: 20
    resources:
        runtime=60 *2,
        mem_mb=mem_allowed,
    shell:
        """
        gappa examine assign --threads {threads} --out-dir {params.outdir} \
            --jplace-path {input.jplace} --taxon-file {input.taxonfile} \
            --ranks-string '{params.ranks_string}' --per-query-results \
            --consensus-thresh {params.consensus_thresh} {params.distribution_ratio} \
            --best-hit --allow-file-overwriting > {log} 2>&1
        """


rule gappa2taxdf:
    """
    Convert gappa output to a taxonomic dataframe
    """
    output:
        "results/{placer}/{ref}/queries/{query}/{heur}/taxonomy.tsv",
    input:
        rules.gappa_assign.output[0],
    log:
        "logs/{placer}/{ref}/queries/{query}/{heur}/gappa2taxdf.log"
    params:
        ranks=lambda wildcards: config["phylogeny"]["ref"][wildcards.ref]["tree_ranks"],
    shell:
        """
        python workflow/scripts/gappa2taxdf.py {input} {output} --ranks {params.ranks} >{log} 2>&1
        """
