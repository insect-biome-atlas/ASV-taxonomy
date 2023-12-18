localrules:
    longest_orfs,

rule run_kaiju:
    input:
        expand("results/kaiju/{ref}/queries/{query}/{query}.classify.tsv",
            ref=config["kaiju"]["ref"].keys(), query=config["kaiju"]["query"].keys())

rule orffinder:
    """
    Translate custom Metabuli reference database to aa for use with kaiju
    """
    output:
        "results/kaiju/{ref}/orffinder/orfs.faa",
    input:
        rules.sintax_fasta_to_gtdbfmt.output.seqs
    log:
        "logs/kaiju/orffinder.{ref}.log",
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
        faa="results/kaiju/{ref}/orffinder/longest.faa",
        txt="results/kaiju/{ref}/orffinder/longest.lengths.txt",
    input:
        fa=rules.orffinder.output,
        taxidmap="results/metabuli/{ref}/db/acc2taxid.map"
    log:
        "logs/kaiju/orffinder/{ref}.longest.log",
    envmodules:
        "bioinfo-tools",
        "biopython"
    params:
        src=srcdir("../scripts/longest_orfs.py")
    shell:
        """
        python {params.src} {input.fa} {output.faa} {output.txt} -s {input.taxidmap} > {log} 2>&1
        """
        
rule kaiju_build:
    output:
        fmi="results/kaiju/{ref}/{ref}.fmi",
        bwt=temp("results/kaiju/{ref}/{ref}.bwt"),
        sa=temp("results/kaiju/{ref}/{ref}.sa")
    input:
        rules.longest_orfs.output.faa
    log:
        mkbwt="logs/kaiju/{ref}/kaiju-mkbwt.log",
        mkfmi="logs/kaiju/{ref}/kaiju-mkfmi.log"
    envmodules:
        "bioinfo-tools",
        "kaiju/1.7.2"
    threads: 20
    resources:
        runtime = 60 * 2,
        mem_mb=mem_allowed
    params:
        db="results/kaiju/{ref}/{ref}"
    shell:
        """
        kaiju-mkbwt -n {threads} -a ACDEFGHIKLMNPQRSTVWY -o {params.db} {input} > {log.mkbwt} 2>&1
        kaiju-mkfmi {params.db} > {log.mkfmi} 2>&1
        """

rule kaiju_classify:
    output:
        "results/kaiju/{ref}/queries/{query}/{query}.classify.tsv"
    input:
        qry=lambda wildcards: config["kaiju"]["query"][wildcards.query],
        fmi=rules.kaiju_build.output.fmi,
        nodes="results/metabuli/{ref}/db/taxonomy/nodes.dmp"
    log:
        "logs/kaiju/kaiju_classify.{ref}.{query}.log"
    params:
        settings=config["kaiju"]["settings"]
    resources:
        runtime = 60 * 10,
        mem_mb=mem_allowed
    envmodules:
        "bioinfo-tools",
        "kaiju/1.7.2"
    threads: 8
    shell:
        """
        kaiju -t {input.nodes} -f {input.fmi} -i {input.qry} -z {threads} -p {params.settings} -o {output} >{log} 2>&1
        """
