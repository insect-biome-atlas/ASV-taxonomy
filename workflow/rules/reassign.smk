def update_taxonomy_input(wildcards):
    return f"results/epa-ng/{wildcards.phyloref}/queries/{wildcards.query}.reassigned/{wildcards.heur}/taxonomy.tsv"

localrules:
    update_taxonomy

rule update_taxonomy:
    output:
        "results/reassign/{ref}/queries/{query}/reassign/{phyloref}/{heur}/taxonomy.tsv"
    input:
        base=rules.parse_sintax.output[0],
        update=update_taxonomy_input,
    params:
        agree_rank = lambda wildcards: config["phylogeny"]["ref"][wildcards.phyloref]["reassign_sintax"]["rank"],
        update_ranks = lambda wildcards: config["phylogeny"]["ref"][wildcards.phyloref]["reassign_sintax"]["update_ranks"],
        downstream_ranks = lambda wildcards: config["phylogeny"]["ref"][wildcards.phyloref]["reassign_sintax"]["downstream_ranks"],
    shell:
        """
        python workflow/scripts/update_taxonomy.py -b {input.base} -u {input.update} \
            -a {params.agree_rank} -U {params.update_ranks} -d {params.downstream_ranks} -o {output}
        """