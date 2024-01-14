import os
import sys

localrules:
    sample_keep_species_in_db,
    sample_keep_species_remove_identical

cases = ["case1-keep-species-in-db", "case2-keep-species-remove-identical", "case3-remove-species-keep-genus",
                     "case4-remove-genus-keep-family", "case5-remove-family"]

# populate config
if "benchmark" in config.keys():
    for db in config["benchmark"].keys():
        ranks = config["benchmark"][db]["ranks"]
        for case in cases:
            short = f"{db}.{case.split("-")[0]}"
            try:
                config["sintax"]["query"][short] = f"benchmark/{db}/{case}/test.fasta"
            except KeyError:
                config["sintax"] = {"query": {short: f"benchmark/{db}/{case}/test.fasta"}}
            try:
                config["sintax"]["ref"][short] = {"fasta": f"benchmark/{db}/{case}/train.fasta"}
            except KeyError:
                config["sintax"] = {"ref": {short: {"fasta": f"benchmark/{db}/{case}/train.fasta"}}}
            config["sintax"]["ref"][short]["ranks"] = ranks    

rule create_testdata:
    input:
        expand("benchmark/{db}/{case}/test.fasta", 
            db=config["benchmark"].keys(), 
            case=cases),

def db_fasta(wildcards):
    try:
        if config["benchmark"][wildcards.db]["extract_subseq"]:            
            return f"benchmark/{wildcards.db}/subseqs.fasta"
    except KeyError:
        pass
    return config["benchmark"][wildcards.db]["fasta"]

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
        fasta = db_fasta,
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
        fasta = db_fasta,
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
        fasta = db_fasta,
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
        fasta = db_fasta,
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
        fasta = db_fasta,
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

## SINTAX EVAL ##
def get_all_sintax_eval(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            input.append(f"results/sintax/{short}/queries/{short}/sintax.eval.tsv")
    return input

rule eval_all_sintax:
    input:
        get_all_sintax_eval

rule evaluate_sintax:
    """
    Evaluate the performance of the classifiers
    """
    output:
        "results/sintax/{ref}/queries/{query}/sintax.eval.tsv"
    input:
        sintax=rules.parse_sintax.output,
        tsv=lambda wildcards: config["sintax"]["query"][wildcards.query].replace(".fasta", ".tsv"),
    log:
        "results/sintax/{ref}/queries/{query}/sintax.eval.log"
    shell:
        """
        python workflow/scripts/evaluate_classifier.py --classifier {input.sintax} --taxfile {input.tsv} -o {output} >{log} 2>&1
        """