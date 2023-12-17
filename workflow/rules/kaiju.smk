rule orffinder:
    """
    Translate custom Metabuli reference database to aa for use with kaiju
    """
    output:
        "results/orffinder/{ref}.orfs.faa",
    input:
        rules.sintax_fasta_to_gtdbfmt.output.seqs
    log:
        "results/orffinder/{ref}.log",
    envmodules:
        "bioinfo-tools",
        "ORFfinder/0.4.3",
    resources:
        runtime = 60 * 4,
        mem_mb = mem_allowed,
    shell:
        """
        ORFfinder -in {input} -ml 30 -g 5 -s 2 -n true -strand plus -out {output} -outfmt 0 > {log} 2>&1
        """


rule longest_orfs:
    output:
        faa="results/orffinder/{ref}.longest.faa",
        txt="results/orffinder/{ref}.longest.lengths.txt",
    input:
        rules.orffinder.output,
    log:
        "logs/orffinder/{ref}.longest.log",
    run:
        seqs = {}
        from Bio.SeqIO import parse
        import re
        for record in parse(input[0], "fasta"):
            seq_id = re.sub("ORF\d+_(\w+).+", r"\g<1>", (record.id).split("|")[1])
            l = len(record.seq)
            if seq_id not in seqs.keys():
                seqs[seq_id] = {"seq": str(record.seq), "id": record.id}
            else:
                if l > len(seqs[seq_id]["seq"]):
                    seqs[seq_id] = {"seq": str(record.seq), "id": record.id}
        with open(output.faa, "w") as fhout, open(output.txt, "w") as fhout_txt:
            for seq_id, d in seqs.items():
                fhout.write(f">{seq_id} {d['id']}\n{d['seq']}\n")
                fhout_txt.write(f"{seq_id}\t{len(d['seq'])}\n")
