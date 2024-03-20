## Phylogenetic datasets ##

localrules:
    sample_keep_species_remove_identical_phylo,
    sample_keep_genus_remove_species_phylo,
    sample_keep_family_remove_genus_phylo,
    sample_keep_order_remove_family_phylo,
    evaluate_phylo,
    eval_all_epang,
    eval_all_pplacer


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
        expand("results/epa-ng/{ref}/queries/{query}/raxml-ng/{heur}/eval.tsv", 
        ref=config["phylogeny"]["ref"].keys(), query=config["phylogeny"]["query"].keys(), heur=config["phylogeny"]["epa-ng"]["heuristics"],
        ),

rule eval_all_pplacer:
    input:
        expand("results/pplacer/{ref}/queries/{query}/{phylotool}/{heur}/eval.tsv", 
        ref=config["phylogeny"]["ref"].keys(), query=config["phylogeny"]["query"].keys(), heur=config["phylogeny"]["pplacer"]["heuristics"],
        phylotool=config["phylogeny"]["pplacer"]["phylo-tools"]
        ),

rule evaluate_phylo:
    output:
        "results/{placer}/{ref}/queries/{query}/{phylotool}/{heur}/eval.tsv"
    input:
        res="results/{placer}/{ref}/queries/{query}/{phylotool}/{heur}/taxonomy.tsv",
        tax=lambda wildcards: config["phylogeny"]["query"][wildcards.query].replace(".fasta", ".tsv"),
    log:
        "logs/{placer}/{ref}/queries/{query}/{phylotool}/{heur}/eval.log"
    params:
        classifier_str = lambda wildcards: wildcards.placer+"."+wildcards.phylotool+"."+wildcards.heur,
    shell:
        """
        python workflow/scripts/evaluate_classifier.py {input.res} --taxonomy {input.tax} --classifier {params.classifier_str} --output {output} >{log} 2>&1
        """