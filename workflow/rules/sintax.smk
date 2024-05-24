localrules:
    parse_sintax,
    run_sintax,
    collate_sintax,
    extract_ASVs

rule run_sintax:
    input: 
        expand("results/sintax/{ref}/queries/{query}/sintax.tsv",
            ref=config["sintax"]["ref"].keys(), query=config["sintax"]["query"].keys())

splits=[f'split{x:03d}' for x in list(range(1,1001))]

rule split_sintax_input:
    """
    Splits the sintax fasta file into 1000 chunks
    """
    output:
        temp(expand("results/sintax/{{ref}}/queries/{{query}}/splits/{split}.fasta", split=splits))
    input:
        qry=lambda wildcards: config["sintax"]["query"][wildcards.query]
    log:
        "logs/sintax/sintax.{ref}.{query}.split.log"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        splits=len(splits),
    resources:
        runtime = 60,
    threads: 2
    shell:
        """
        cat {input.qry} | seqkit split2 -O {params.outdir} -j {threads} -p {params.splits} --by-part-prefix split >{log} 2>&1
        """

rule sintax:
    """
    Runs sintax on one split of the query fasta file
    """
    output:
        temp("results/sintax/{ref}/queries/{query}/splits/{split}.tab")
    input:
        db=lambda wildcards: config["sintax"]["ref"][wildcards.ref]["fasta"],
        qry="results/sintax/{ref}/queries/{query}/splits/{split}.fasta"
    log:
        "logs/sintax/splits/sintax.{ref}.{query}.{split}.log"
    params:
        seed=config["sintax"]["randseed"],
        cutoff=config["sintax"]["cutoff"]
    conda: "../envs/vsearch.yml"
    resources:
        runtime = 30,
        constraint="mem256GB"
    threads: 20
    shell:
        """
        vsearch --sintax {input.qry} --sintax_cutoff {params.cutoff} --randseed {params.seed} --db {input.db} --tabbedout {output} --threads 1 >{log} 2>&1
        """

rule collate_sintax:
    """
    Concatenates the sintax output files into a single file
    """
    output:
        "results/sintax/{ref}/queries/{query}/sintax.tab"
    input:
        expand("results/sintax/{{ref}}/queries/{{query}}/splits/{split}.tab", split=splits),
    shell:
        """
        cat {input} > {output}
        """

rule parse_sintax:
    """
    Parses the sintax output file into a tsv file
    """
    output:
        "results/sintax/{ref}/queries/{query}/sintax.tsv"
    input:
        rules.collate_sintax.output
    log:
        "logs/sintax/parse_sintax.{ref}.{query}.log"
    params:
        src="workflow/scripts/sintax_tsv.py",
        ranks=lambda wildcards: config["sintax"]["ref"][wildcards.ref]["ranks"]
    shell:
        """
        python {params.src} -i {input} -o {output} -r {params.ranks} > {log} 2>&1
        """

rule extract_ASVs:
    """
    Extract ASVs matching reassign_rank == reassign_taxa but unclassified at child rank.
    For example, ASVs classified as reassign_taxa == 'Insecta' at reassign_rank == 'Class' but unclassified at 'Order'.
    """
    output:
        tsv="results/sintax/{ref}/queries/{query}/reassign/{rank}/unclassified.{taxa}.tsv",
        fasta="results/sintax/{ref}/queries/{query}/reassign/{rank}/unclassified.{taxa}.fasta"
    input:
        tsv=rules.parse_sintax.output[0],
        qry=lambda wildcards: config["sintax"]["query"][wildcards.query]
    params:
        reassign_rank = lambda wildcards: wildcards.rank,
        reassign_taxa = lambda wildcards: wildcards.taxa,
    run:
        reassign_taxa = params.reassign_taxa
        reassign_rank = params.reassign_rank
        from Bio.SeqIO import parse, write as write_fasta
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t", index_col=0)
        cols = df.columns.tolist()
        child_rank = cols[cols.index(reassign_rank)+1]
        taxdf = df.loc[df[child_rank]==f"unclassified.{reassign_taxa}"]
        seqs = []
        for record in parse(input.qry, "fasta"):
            if record.id in taxdf.index:
                seqs.append(record)
        with open(output.fasta, "w") as f:
            write_fasta(seqs, f, "fasta")
        taxdf.to_csv(output.tsv, sep="\t")