## Phylogenetic datasets ##

localrules:
    sample_keep_species_remove_identical_phylo,
    sample_keep_genus_remove_species_phylo,
    sample_keep_family_remove_genus_phylo,
    sample_keep_order_remove_family_phylo,
    evaluate_epang,
    eval_all_epang


phylo_cases = ["case1-sample_keep_species_remove_identical_phylo",
               "case2-sample_keep_genus_remove_species_phylo",
               "case3-sample_keep_family_remove_genus_phylo",
               "case4-sample_keep_order_remove_family_phylo"]

rule create_phylo_testdata:
    input:
        expand("benchmark_phylo/{db}/{case}/test.fasta", 
            db=config["benchmark_phylo"].keys(), 
            case=phylo_cases),

rule sample_keep_species_remove_identical_phylo:
    """
    Easy positive case: coidb sequences for which species, but not identical sequences, exist in tree [1]
    """
    output:
        test_fasta = "benchmark_phylo/{db}/case1-sample_keep_species_remove_identical_phylo/test.fasta",
        test_tsv = "benchmark_phylo/{db}/case1-sample_keep_species_remove_identical_phylo/test.tsv",
    input:
        tree_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_fasta"],
        db_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["fasta"],
        tree_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_taxonomy"],
        db_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["taxonomy"],
    log:
        "benchmark_phylo/{db}/case1-sample_keep_species_remove_identical_phylo/log.txt",
    params:
        k = 100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv),
        ranks = lambda wildcards: config["benchmark_phylo"][wildcards.db]["ranks"],
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
    threads: 2
    shell:
        """
        python workflow/scripts/create_testdata_phylo.py \
            --input_db_fasta {input.db_fasta} \
            --input_tree_fasta {input.tree_fasta} \
            --input_tree_taxfile {input.tree_taxonomy} \
            --input_db_taxfile {input.db_taxonomy} \
            --output_dir {params.outdir} \
            -k {params.k} \
            -s {params.seed} \
            --threads {threads} --ranks {params.ranks} \
            --case 1 > {log} 2>&1
        """

rule sample_keep_genus_remove_species_phylo:
    """
    coidb sequences for which genus, but not species, exist in tree [2]
    """
    output:
        test_fasta = "benchmark_phylo/{db}/case2-sample_keep_genus_remove_species_phylo/test.fasta",
        test_tsv = "benchmark_phylo/{db}/case2-sample_keep_genus_remove_species_phylo/test.tsv",
    input:
        tree_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_fasta"],
        db_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["fasta"],
        tree_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_taxonomy"],
        db_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["taxonomy"],
    log:
        "benchmark_phylo/{db}/case2-sample_keep_genus_remove_species_phylo/log.txt",
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv),
        ranks=lambda wildcards: config["benchmark_phylo"][wildcards.db]["ranks"],
    threads: 1
    shell:
        """
        python workflow/scripts/create_testdata_phylo.py \
            --input_db_fasta {input.db_fasta} \
            --input_tree_fasta {input.tree_fasta} \
            --input_tree_taxfile {input.tree_taxonomy} \
            --input_db_taxfile {input.db_taxonomy} \
            --output_dir {params.outdir} \
            -k {params.k} \
            -s {params.seed} \
            --threads {threads} \
            --case 2 --ranks {params.ranks} > {log} 2>&1
        """

rule sample_keep_family_remove_genus_phylo:
    """
    coidb sequences for which family, but not genus, exist in tree [3]
    """
    output:
        test_fasta = "benchmark_phylo/{db}/case3-sample_keep_family_remove_genus_phylo/test.fasta",
        test_tsv = "benchmark_phylo/{db}/case3-sample_keep_family_remove_genus_phylo/test.tsv",
    input:
        tree_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_fasta"],
        db_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["fasta"],
        tree_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_taxonomy"],
        db_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["taxonomy"],
    log:
        "benchmark_phylo/{db}/case3-sample_keep_family_remove_genus_phylo/log.txt",
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv),
        ranks=lambda wildcards: config["benchmark_phylo"][wildcards.db]["ranks"],
    threads: 1
    shell:
        """
        python workflow/scripts/create_testdata_phylo.py \
            --input_db_fasta {input.db_fasta} \
            --input_tree_fasta {input.tree_fasta} \
            --input_tree_taxfile {input.tree_taxonomy} \
            --input_db_taxfile {input.db_taxonomy} \
            --output_dir {params.outdir} \
            -k {params.k} \
            -s {params.seed} \
            --threads {threads} \
            --case 3 --ranks {params.ranks} > {log} 2>&1
        """

rule sample_keep_order_remove_family_phylo:
    """
    coidb sequences for which order, but not family, exist in tree [4]
    """
    output:
        test_fasta = "benchmark_phylo/{db}/case4-sample_keep_order_remove_family_phylo/test.fasta",
        test_tsv = "benchmark_phylo/{db}/case4-sample_keep_order_remove_family_phylo/test.tsv",
    input:
        tree_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_fasta"],
        db_fasta = lambda wildcards: config["benchmark_phylo"][wildcards.db]["fasta"],
        tree_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["tree_taxonomy"],
        db_taxonomy = lambda wildcards: config["benchmark_phylo"][wildcards.db]["taxonomy"],
    log:
        "benchmark_phylo/{db}/case4-sample_keep_order_remove_family_phylo/log.txt",
    params:
        k=100,
        seed=42,
        outdir=lambda wildcards, output: os.path.dirname(output.test_tsv),
        ranks=lambda wildcards: config["benchmark_phylo"][wildcards.db]["ranks"],
    threads: 1
    shell:
        """
        python workflow/scripts/create_testdata_phylo.py \
            --input_db_fasta {input.db_fasta} \
            --input_tree_fasta {input.tree_fasta} \
            --input_tree_taxfile {input.tree_taxonomy} \
            --input_db_taxfile {input.db_taxonomy} \
            --output_dir {params.outdir} \
            -k {params.k} \
            -s {params.seed} \
            --threads {threads} \
            --case 4 --ranks {params.ranks} > {log} 2>&1
        """

rule eval_all_epang:
    input:
        expand("results/epa-ng/{ref}/queries/{query}/{query}.eval.tsv", 
        ref=config["epa-ng"]["ref"].keys(), query=config["epa-ng"]["query"].keys()
        ),

rule evaluate_epang:
    output:
        "results/epa-ng/{ref}/queries/{query}/{query}.eval.tsv"
    input:
        res="results/epa-ng/{ref}/queries/{query}.taxonomy.tsv",
        tax=lambda wildcards: config["epa-ng"]["query"][wildcards.query].replace(".fasta", ".tsv"),
    log:
        "results/epa-ng/{ref}/queries/{query}/epa-ng.eval.log"
    params:
        classifier_str = "epa-ng",
    shell:
        """
        python workflow/scripts/evaluate_classifier.py {input.res} --taxonomy {input.tax} --classifier {params.classifier_str} --output {output} >{log} 2>&1
        """