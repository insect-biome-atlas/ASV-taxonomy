import os

localrules:
    sample_keep_species_in_db,
    sample_keep_species_remove_identical

rule create_testdata:
    input:
        expand("benchmark/{db}/{case}/test.fasta", 
            db=config["benchmark"].keys(), 
            case=["case1-keep-species-in-db", "case2-keep-species-remove-identical", "case3-remove-species-keep-genus",
                  "case4-remove-genus-keep-family", "case5-remove-family"]),

rule mafft_align_db:
    """
    Align the database sequences
    """
    output:
        "benchmark/{db}/mafft-aligned.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
    log:
        "benchmark/{db}/mafft-align_db.log"
    threads: 20
    envmodules:
        "bioinfo-tools",
        "MAFFT/7.407"
    resources:
        runtime = 60 * 24 * 10, # 10 days
        mem_mb = mem_allowed,
    shell:
        """
        mafft --thread {threads} {input.fasta} >{output} 2>{log}
        """

rule sample_keep_species_in_db:
    """
    For this case, the train fasta is the same as the original database file, so use symlink
    """
    output:
        train_tsv="benchmark/{db}/case1-keep-species-in-db/train.tsv",
        train_fasta="benchmark/{db}/case1-keep-species-in-db/train.fasta",
        test_tsv="benchmark/{db}/case1-keep-species-in-db/test.tsv",
        test_fasta="benchmark/{db}/case1-keep-species-in-db/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    log:
        "benchmark/{db}/case1-keep-species-in-db/log.txt"
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv)
    shell:
        """
        python workflow/scripts/create_testdata.py --input_fasta {input.fasta} \
            --input_taxfile {input.tax} \
            --output_dir {params.outdir} \
            -k {params.k} -s {params.seed} --case 1 >{log} 2>&1
        """

rule sample_keep_species_remove_identical:
    """
    Hard positive case: keep species in db but remove identical sequences.
    Sample 100 species, for each species dereplicate to unique sequences, then sample 1 sequence 
    from the dereplicated sequences and remove any identical sequences from the db.
    """
    output:
        train_tsv="benchmark/{db}/case2-keep-species-remove-identical/train.tsv",
        train_fasta="benchmark/{db}/case2-keep-species-remove-identical/train.fasta",
        test_tsv="benchmark/{db}/case2-keep-species-remove-identical/test.tsv",
        test_fasta="benchmark/{db}/case2-keep-species-remove-identical/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    log:
        "benchmark/{db}/case2-keep-species-remove-identical/log.txt"
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv)
    shell:
        """
        python workflow/scripts/create_testdata.py --input_fasta {input.fasta} \
            --input_taxfile {input.tax} \
            --output_dir {params.outdir} \
            -k {params.k} -s {params.seed} --case 2 >{log} 2>&1
        """
        
rule sample_keep_genus:
    """
    Negative case: remove species in db, but make sure at least 1 species for the genus is present
    """
    output:
        train_tsv="benchmark/{db}/case3-remove-species-keep-genus/train.tsv",
        train_fasta="benchmark/{db}/case3-remove-species-keep-genus/train.fasta",
        test_tsv="benchmark/{db}/case3-remove-species-keep-genus/test.tsv",
        test_fasta="benchmark/{db}/case3-remove-species-keep-genus/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    log:
        "benchmark/{db}/case3-remove-species-keep-genus/log.txt"
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv)
    shell:
        """
        python workflow/scripts/create_testdata.py --input_fasta {input.fasta} \
            --input_taxfile {input.tax} \
            --output_dir {params.outdir} \
            -k {params.k} -s {params.seed} --case 3 >{log} 2>&1
        """

rule sample_keep_family:
    """
    Negative case: remove genus in db, but make sure at least 1 genus for the family is present.
    """
    output:
        train_tsv="benchmark/{db}/case4-remove-genus-keep-family/train.tsv",
        train_fasta="benchmark/{db}/case4-remove-genus-keep-family/train.fasta",
        test_tsv="benchmark/{db}/case4-remove-genus-keep-family/test.tsv",
        test_fasta="benchmark/{db}/case4-remove-genus-keep-family/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    log:
        "benchmark/{db}/case4-remove-genus-keep-family/log.txt"
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv)
    shell:
        """
        python workflow/scripts/create_testdata.py --input_fasta {input.fasta} \
            --input_taxfile {input.tax} \
            --output_dir {params.outdir} \
            -k {params.k} -s {params.seed} --case 4 >{log} 2>&1
        """

rule sample_remove_family:
    """
    Negative case: remove all sequences for the family
    """
    output:
        train_tsv="benchmark/{db}/case5-remove-family/train.tsv",
        train_fasta="benchmark/{db}/case5-remove-family/train.fasta",
        test_tsv="benchmark/{db}/case5-remove-family/test.tsv",
        test_fasta="benchmark/{db}/case5-remove-family/test.fasta",
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        tax = lambda wildcards: config["benchmark"][wildcards.db]["taxonomy"],
    log:
        "benchmark/{db}/case5-remove-family/log.txt"
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv)
    shell:
        """
        python workflow/scripts/create_testdata.py --input_fasta {input.fasta} \
            --input_taxfile {input.tax} \
            --output_dir {params.outdir} \
            -k {params.k} -s {params.seed} --case 5 >{log} 2>&1
        """