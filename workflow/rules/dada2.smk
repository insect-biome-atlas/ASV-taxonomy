wildcard_constraints:
    boot="[0-9]+",

localrules:
    sintax2dada2,

rule run_dada2:
    input:
        expand("results/dada2/{ref}/queries/{query}/dada2.minBoot{boot}.tsv", 
                ref=config["dada2"]["ref"], query=config["dada2"]["query"], boot = config["dada2"]["minBoot"])
        
rule sintax2dada2:
    """
    Converts sintax formatted fasta to assignTaxonomy format
    """
    input:
        fasta = lambda wildcards: config["dada2"]["ref"][wildcards.ref]["fasta"],
    output:
        fasta = "results/dada2/{ref}/assignTaxonomy.fasta",
    params:
        ranks = lambda wildcards: config["dada2"]["ref"][wildcards.ref]["ranks"]
    run:
        from Bio.SeqIO import parse
        import re
        ranks = [x.lower() for x in params.ranks]
        with open(input.fasta, 'r') as fhin, open(output.fasta, 'w') as fhout:
            for record in parse(fhin, "fasta"):
                desc = ";".join((record.description).split("=")[-1].split(",")[0:len(ranks)])
                desc = re.sub('[a-z]:', '', desc+";")
                fhout.write(f">{desc}\n{str(record.seq)}\n")

def get_dada2_ref(wildcards):
    if config["dada2"]["ref"][wildcards.ref]["format"] == "sintax":
        return f"results/dada2/{wildcards.ref}/assignTaxonomy.fasta"
    else:
        return config["dada2"]["ref"][wildcards.ref]["fasta"]

rule dada2:
    output:
        "results/dada2/{ref}/queries/{query}/dada2.minBoot{boot}.tsv"
    input:
        ref=get_dada2_ref,
        qry=lambda wildcards: config["dada2"]["query"][wildcards.query]
    log:
        "results/dada2/{ref}/queries/{query}/dada2.minBoot{boot}.log"
    params:
        taxLevels = lambda wildcards: config["dada2"]["ref"][wildcards.ref]["ranks"],
        seed = 100,
    threads: 10
    conda:
        "../envs/dada2.yml"
    resources:
        mem_mb=mem_allowed,
        runtime = 60 * 10
    #container:
    #    "docker://quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0"
    script:
        "../scripts/dada2.R"