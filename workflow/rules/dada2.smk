rule sintax2dada2:
    """
    Converts sintax formatted fasta to assignTaxonomy format
    Truncates taxLevels to 7
    """
    input:
        fasta = lambda wildcards: config["dada2"]["ref"][wildcards.ref]["fasta"],
    output:
        fasta = "results/dada2/{ref}/assignTaxonomy.fasta",
    run:
        from Bio.SeqIO import parse
        import re
        with open(input.fasta, 'r') as fhin, open(output.fasta, 'w') as fhout:
            for record in parse(fhin, "fasta"):
                desc = ";".join((record.description).split("=")[-1].split(",")[0:7])
                desc = re.sub('[a-z]:', '', desc)+";"
                fhout.write(f">{desc}\n{str(record.seq)}\n")

def get_dada2_ref(wildcards):
    if config["dada2"]["ref"][wildcards.ref]["fmt"] == "sintax":
        return f"results/dada2/{wildcards.ref}/assignTaxonomy.fasta"
    else:
        return config["dada2"]["ref"][wildcards.ref]["fasta"]

rule dada2:
    output:
        "results/dada2/{ref}/queries/{query}/dada2.tsv"
    input:
        ref=get_dada2_ref,
        qry=lambda wildcards: config["dada2"]["query"][wildcards.query]
    params:
        taxLevels = lambda wildcards: config["dada2"]["ref"][wildcards.ref]["ranks"],
        seed = 42,
    threads: 4
    conda:
        "../envs/dada2.yml"
    script:
        "../scripts/dada2.R"

    