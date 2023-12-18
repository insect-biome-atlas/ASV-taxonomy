localrules:
    parse_sintax,
    run_sintax

rule run_sintax:
    input: 
        expand("results/sintax/{ref}/queries/{query}/sintax.tsv",
            ref=config["sintax"]["ref"].keys(), query=config["sintax"]["query"].keys())

rule sintax:
    output:
        "results/sintax/{ref}/queries/{query}/sintax.tab"
    input:
        db=lambda wildcards: config["sintax"]["ref"][wildcards.ref],
        qry=lambda wildcards: config["sintax"]["query"][wildcards.query]
    log:
        "logs/sintax/sintax.{ref}.{query}.log"
    params:
        seed=config["sintax"]["randseed"],
        cutoff=config["sintax"]["cutoff"]
    conda: "../envs/vsearch.yml"
    resources:
        runtime = 60 * 24
    threads: 4
    shell:
        """
        vsearch --sintax {input.qry} --sintax_cutoff {params.cutoff} --randseed {params.seed} --db {input.db} --tabbedout {output} --threads 1 >{log} 2>&1
        """

rule parse_sintax:
    output:
        "results/sintax/{ref}/queries/{query}/sintax.tsv"
    input:
        rules.sintax.output
    log:
        "logs/sintax/parse_sintax.{ref}.{query}.log"
    params:
        src=srcdir("../scripts/sintax_tsv.py")
    shell:
        """
        python {params.src} -i {input} -o {output} > {log} 2>&1
        """