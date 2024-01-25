localrules:
    sintax2qiime_input,

rule sintax2qiime_input:
    """
    Output should look like this:
    seqid1	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
    seqid2	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__

    fasta output should be:
    >seqid1
    ATGCGGGCTAGAGTAGCGAT...
    """
    input:
        tsv=lambda wildcards: config["qiime2"]["ref"][wildcards.ref]["taxonomy"],
        fasta=lambda wildcards: config["qiime2"]["ref"][wildcards.ref]["fasta"],
    output:
        tsv="results/qiime2/{ref}/taxonomy.tsv",
        fasta="results/qiime2/{ref}/seqs.fasta",
    run:
        import pandas as pd
        from Bio.SeqIO import parse
        df = pd.read_csv(input.tsv, sep="\t", index_col=0)
        seqs = {}
        for record in parse(input.fasta, "fasta"):
            seqid = (record.id).split(";")[0]
            seqs[seqid] = record
        common_seqs = list(set(seqs.keys()).intersection(list(df.index)))
        df = df.loc[common_seqs]
        df.to_csv(output.tsv, sep="\t")
        with open(output.fasta, 'r') as fhout:
            for seqid, record in seqs.items():
            fhout.write(f">{seqid}\n{str(record.seq)}\n")

        