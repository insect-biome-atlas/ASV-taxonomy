localrules:
    nexus2newick,
    extract_ref_taxonomy,
    nexus2fasta,
    split_epang_input,
    split_aln,
    gappa2taxdf,
    collate_gappa_taxdf,
    write_config,
    write_software,

wildcard_constraints:
    heur="baseball-heur|dyn-heur|no-heur",
    placer="epa-ng",

## target rules

def epa_ng_input(wildcards):
    input = []
    for ref in config["phylogeny"]["ref"].keys():
        if "query" in config["phylogeny"].keys():
            for query in config["phylogeny"]["query"].keys():
                for heur in config["phylogeny"]["epa-ng"]["heuristics"]:
                    input.append(f"results/epa-ng/{ref}/queries/{query}/{heur}/taxonomy.tsv")
        if "reassign_sintax" in config["phylogeny"]["ref"][ref].keys():
            rank = config["phylogeny"]["ref"][ref]["reassign_sintax"]["rank"]
            taxa = config["phylogeny"]["ref"][ref]["reassign_sintax"]["taxa"]
            for sintax_ref in config["sintax"]["ref"].keys():
                for sintax_query in config["sintax"]["query"].keys():
                    new_query = f"{sintax_query}.reassigned"
                    for heur in config["phylogeny"]["epa-ng"]["heuristics"]:
                        input.append(f"results/epa-ng/{ref}/queries/{new_query}/{heur}/taxonomy.tsv")
                        input.append(f"results/reassign/{sintax_ref}/queries/{sintax_query}/reassign/{ref}/{heur}/taxonomy.tsv")
    return input

rule run_epa_ng:
    input:
        epa_ng_input,

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

checkpoint split_epang_input:
    """
    Splits the sintax fasta file into chunks with fixed size
    """
    output:
        directory("results/phylogeny/{ref}/queries/{query}/splits")
    input:
        qry=lambda wildcards: config["phylogeny"]["query"][wildcards.query],
    log:
        "logs/phylogeny/{ref}.{query}.split.log"
    params:
        outdir=lambda wildcards, output: output[0],
        size=500
    resources:
        runtime = 60,
    threads: 1
    shell:
        """
        cat {input.qry} | seqkit split2 -O {params.outdir} -j {threads} -s {params.size} --by-part-prefix split >{log} 2>&1
        """

rule hmm_align:
    output:
        temp("results/phylogeny/{ref}/hmmalign/{query}/splits/{split}/{split}.fasta"),
    input:
        hmm=rules.hmm_build.output,
        qry="results/phylogeny/{ref}/queries/{query}/splits/stdin.part_{split}.fasta",
        ref_msa=ref_msa
    log:
        "logs/phylogeny/{ref}/hmmalign.{query}.{split}.log",
    resources:
        runtime=60,
        mem_mb=mem_allowed,
    params:
        tmpdir=lambda wildcards: f"$TMPDIR/{wildcards.ref}.{wildcards.query}.{wildcards.split}.hmm_align",
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
        ref_msa=temp("results/phylogeny/{ref}/hmmalign/{query}/splits/{split}/reference.fasta"),
        qry_msa=temp("results/phylogeny/{ref}/hmmalign/{query}/splits/{split}/query.fasta"),
    input:
        ref_msa=ref_msa,
        msa=rules.hmm_align.output,
    log:
        "logs/phylogeny/{ref}/split_aln.{query}.{split}.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.ref_msa),
    conda: "../envs/epang.yml"
    shell:
        """
        epa-ng --redo --split {input.ref_msa} {input.msa} --outdir {params.outdir} > {log} 2>&1
        """

rule raxml_evaluate:
    output:
        temp("results/phylogeny/{ref}/raxml-ng/{query}/splits/{split}/info.raxml.bestModel"),
    input:
        tree=ref_tree,
        msa=rules.split_aln.output.ref_msa,
    log:
        "logs/phylogeny/{ref}/raxml-ng.{query}.{split}.log",
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
        temp("results/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/epa-ng_result.jplace"),
    input:
        qry=rules.split_aln.output.qry_msa,
        ref_msa=rules.split_aln.output.ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output[0],
    log:
        "logs/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/epa-ng.log",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        heur=get_heuristic,
    conda: "../envs/epang.yml"
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
    """
    output:
        temp("results/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/per_query.tsv"),
        temp("results/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/profile.tsv"),
        temp("results/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/labelled_tree.newick"),
    input:
        jplace=rules.epa_ng.output,
        taxonfile=ref_taxonomy,
    log:
        "logs/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/gappa_assign.log",
    params:
        ranks_string=lambda wildcards: "|".join(config["phylogeny"]["ref"][wildcards.ref]["tree_ranks"]),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        consensus_thresh=config["phylogeny"]["gappa"]["consensus_thresh"],
        distribution_ratio=get_dist_ratio(config),
    conda:
        "../envs/gappa.yml"
    threads: 20
    resources:
        runtime=30,
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
        temp("results/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/taxonomy.tsv"),
    input:
        rules.gappa_assign.output[0],
    log:
        "logs/epa-ng/{ref}/queries/{query}/splits/{heur}/gappa2taxdf.{split}.log"
    params:
        ranks=lambda wildcards: config["phylogeny"]["ref"][wildcards.ref]["tree_ranks"],
    shell:
        """
        python workflow/scripts/gappa2taxdf.py {input} {output} --ranks {params.ranks} >{log} 2>&1
        """

def aggregate_gappa(wildcards):
    checkpoint_output = checkpoints.split_epang_input.get(**wildcards).output[0]
    return expand("results/epa-ng/{ref}/queries/{query}/splits/{split}/{heur}/taxonomy.tsv",
                    ref=wildcards.ref, query=wildcards.query, heur=wildcards.heur, 
                    split=glob_wildcards(os.path.join(checkpoint_output, "stdin.part_{split}.fasta")).split)


rule collate_gappa_taxdf:
    output:
        "results/epa-ng/{ref}/queries/{query}/{heur}/taxonomy.tsv"
    input:
        aggregate_gappa
    run:
        with open(output[0], "w") as out:
            for i, f in enumerate(input):
                with open(f, 'r') as infile:
                    for j, line in enumerate(infile):
                        if j == 0 and i == 0:
                            out.write(line)
                        elif j > 0:
                            out.write(line)