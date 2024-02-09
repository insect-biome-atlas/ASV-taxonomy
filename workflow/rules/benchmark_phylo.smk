## Phylogenetic datasets ##

rule sample_keep_species_remove_identical_phylo:
    """
    Easy positive case: coidb sequences for which species, but not identical sequences, exist in tree [1]
    """
    output:
        test_fasta = "results/benchmark_phylo/case1-sample_keep_species_remove_identical_phylo/test.fasta",
        test_tsv = "results/benchmark_phylo/case1-sample_keep_species_remove_identical_phylo/test.tsv",
    input:
        tree_fasta = config["benchmark_phylo"][wildcards.db]["tree_fasta"],
        db_fasta = config["benchmark_phylo"][wildcards.db]["fasta"],
        tree_taxonomy = config["benchmark_phylo"][wildcards.db]["tree_taxonomy"],
        db_taxonomy = config["benchmark_phylo"][wildcards.db]["taxonomy"],
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.test_fasta),
        k = 100,
    threads: 1
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
    shell:
        """
        python workflow/scripts/create_testdata_phylo.py \
            --input_db_fasta {input.db_fasta} \
            --input_tree_fasta {input.tree_fasta} \
            --input_tree_taxfile {input.tree_taxonomy} \
            --input_db_taxfile {input.db_taxonomy} \
            --output_dir {params.output_dir} \
            -k {params.k} \
            -s {params.seed} \
            --threads {threads} \
            --case 1 > {log} 2>&1
        """

