localrules:
    nexus2newick,
    extract_ref_taxonomy,
    nexus2fasta,
    split_aln,
    gappa2taxdf,
    write_config,
    write_software,

wildcard_constraints:
    heur="baseball-heur|dyn-heur|no-heur"

rule run_epa_ng:
    input:
        expand("results/epa-ng/{ref}/queries/{query}/{heur}/taxonomy.tsv",
            ref = config["epa-ng"]["ref"].keys(), query = config["epa-ng"]["query"].keys(), heur = config["epa-ng"]["heuristics"])

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

def get_heuristic(wildcards):
    if wildcards.heur == "no-heur":
        return "--no-heur"
    elif wildcards.heur == "dyn-heur":
        return "--dyn-heur 0.99999"
    elif wildcards.heur == "baseball-heur":
        return "--baseball-heur"

rule epa_ng:
    output:
        "results/epa-ng/{ref}/queries/{query}/{heur}/epa_result.jplace",
    input:
        qry=rules.split_aln.output.qry_msa,
        ref_msa=rules.split_aln.output.ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output
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
        """

rule taxit_create:
    output:
        "results/pplacer/{ref}/taxit/refpkg/CONTENTS.json"
    input:
        ref_msa=ref_msa,
        ref_tree=ref_tree,
        info=rules.raxml_evaluate.output,
    log:
        "logs/taxit/{ref}/create.log"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        taxit create -l {wildcards.ref} -P {params.outdir} --aln-fasta {input.ref_msa} --tree-stats {input.info} -tree-file {input.ref_tree} 2> {log}
        """

rule pplacer:
    output:
        jplace="results/pplacer/{ref}/queries/{query}/{query}.jplace"
    input:
        repkg=rules.taxit_create.output[0],
        msa=rules.hmm_align.output[0]
    log:
        "logs/pplacer/{ref}/{query}.log"
    shell:
        """
        pplacer -c {input.refpkg} {input.msa} -o {output.jplace} --timing > {log} 2>&1
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
        "results/epa-ng/{ref}/queries/{query}/{heur}/per_query.tsv",
        "results/epa-ng/{ref}/queries/{query}/{heur}/profile.tsv",
        "results/epa-ng/{ref}/queries/{query}/{heur}/labelled_tree.newick",
    input:
        json=rules.epa_ng.output,
        taxonfile=ref_taxonomy,
    log:
        "logs/epa-ng/{ref}/queries/{query}/{heur}/gappa_assign.log",
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
        "results/epa-ng/{ref}/queries/{query}/{heur}/taxonomy.tsv",
    input:
        rules.gappa_assign.output[0],
    log:
        "logs/epa-ng/{ref}/queries/{query}/{heur}/gappa2taxdf.log"
    params:
        ranks=lambda wildcards: config["epa-ng"]["ref"][wildcards.ref]["tree_ranks"],
    shell:
        """
        python workflow/scripts/gappa2taxdf.py {input} {output} --ranks {params.ranks} >{log} 2>&1
        """

## plusplacer-taxtastic

rule fasttree_info:
    output:
        info="results/plusplacer-taxtastic/{ref}/fasttree/info.txt",
        tree="results/plusplacer-taxtastic/{ref}/fasttree/best-tree.nwk"
    input:
        ref_tree = lambda wildcards: config["plusplacer-taxtastic"]["ref"][wildcards.ref]["tree"],
        ref_msa = lambda wildcards: config["plusplacer-taxtastic"]["ref"][wildcards.ref]["msa"],
    log:
        "logs/plusplacer-taxtastic/{ref}/fasttree/fasttree_info.log"
    threads: 1
    resources:
        runtime = 60*24
    conda: "../envs/plusplacer-taxtastic.yml"
    shell:
        """
        FastTree -nosupport -gtr -gamma -nt -log {output.info} -intree {input.ref_tree} < {input.ref_msa} > {output.tree} 2>{log}
        """

rule pplacer_tax_SCAMPP:
    output:
        jplace="results/plusplacer-taxtastic/{ref}/queries/{query}/output.jplace"
    input:
        info=rules.fasttree_info.output.info,
        alignment="results/epa-ng/{ref}/hmmalign/{query}/{ref}.{query}.fasta", # alignment with both queries and references
        ref_msa=ref_msa,
        tree=rules.fasttree_info.output.tree,
    log:
        "logs/plusplacer-taxtastic/{ref}/queries/{query}/pplacer_tax_SCAMPP.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.jplace),
        outname = lambda wildcards, output: os.path.basename(output.jplace).replace(".jplace", ""),
    shadow: "minimal"
    threads: 1
    resources:
        runtime = 60*24
    conda: "../envs/plusplacer-taxtastic.yml"
    shell:
        """
        pplacer-tax-SCAMPP -i {input.info} -d {params.outdir} -r {input.ref_msa} -a {input.alignment} -t {input.tree} -o {params.outname} >{log} 2>&1
        """