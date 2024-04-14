import os
import sys


localrules:
    evaluate_dada2,
    evaluate_sintax,
    evaluate_qiime2,
    eval_all_dada2,
    eval_all_sintax,
    eval_all_qiime2,
    eval_all_qiime2_sklearn,
    eval_all_qiime2_vsearch

cases = ["case1-keep-species-in-db", "case2-keep-species-remove-identical", "case3-remove-species-keep-genus",
                     "case4-remove-genus-keep-family", "case5-remove-family"]

rule create_testdata:
    input:
        expand("benchmark/{db}/{case}/test.fasta", 
            db=config["benchmark"].keys(), 
            case=cases),

def db_fasta(wildcards):
    try:
        if config["benchmark"][wildcards.db]["extract_subseq"]:
            include: "subseq.smk"
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
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
    threads: 2
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
    threads: 2
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
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
    threads: 2
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
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
    threads: 2
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
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
    threads: 2
    resources:
        runtime = 60,
        mem_mb = mem_allowed,
    shell:
        """
        python workflow/scripts/create_testdata.py --input_fasta {input.fasta} \
            --input_taxfile {input.tax} \
            --output_dir {params.outdir} \
            -k {params.k} -s {params.seed} --case 5 >{log} 2>&1
        """

## SINTAX EVAL ##
def get_sintax_benchmark(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            input.append(f"results/sintax/{short}/queries/{short}/sintax.tsv")
    return input

rule benchmark_sintax:
    input:
        get_sintax_benchmark

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
    Evaluate the performance of SINTAX
    """
    output:
        "results/sintax/{ref}/queries/{query}/sintax.eval.tsv"
    input:
        res=rules.parse_sintax.output,
        tax=lambda wildcards: config["sintax"]["query"][wildcards.query].replace(".fasta", ".tsv"),
    log:
        "results/sintax/{ref}/queries/{query}/sintax.eval.log"
    shell:
        """
        python workflow/scripts/evaluate_classifier.py {input.res} --taxonomy {input.tax} --classifier sintax --output {output} >{log} 2>&1
        """

## DADA2 EVAL ##
def get_dada2_benchmark(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            input += expand("results/dada2/{short}/queries/{short}/dada2.minBoot{boot}.tsv",
                        boot=config["dada2"]["minBoot"], short=short)
    return input

rule benchmark_dada2:
    input:
        get_dada2_benchmark

def get_all_dada2_eval(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            input +=expand("results/dada2/{short}/queries/{short}/dada2.minBoot{boot}.eval.tsv",
                        boot=config["dada2"]["minBoot"], short=short)
    return input

rule eval_all_dada2:
    input:
        get_all_dada2_eval

rule evaluate_dada2:
    """
    Evaluate the performance of DADA2
    """
    output:
        "results/dada2/{ref}/queries/{query}/dada2.minBoot{boot}.eval.tsv"
    input:
        res="results/dada2/{ref}/queries/{query}/dada2.minBoot{boot}.tsv",
        tax=lambda wildcards: config["dada2"]["query"][wildcards.query].replace(".fasta", ".tsv"),
    log:
        "results/dada2/{ref}/queries/{query}/dada2.eval.minBoot{boot}.log"
    params:
        classifier_str = lambda wildcards: f"dada2.minBoot{wildcards.boot}"
    shell:
        """
        python workflow/scripts/evaluate_classifier.py {input.res} --taxonomy {input.tax} --classifier {params.classifier_str} --output {output} >{log} 2>&1
        """

## QIIME2 EVAL ##
def get_qiime2_benchmark(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            for classifier in ["sklearn", "vsearch"]:
                input.append(f"results/qiime2/{short}/queries/{short}/taxonomy_{classifier}.tsv")
    return input

rule benchmark_qiime2:
    input:
        get_qiime2_benchmark

def get_all_qiime2_eval(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            for classifier in ["sklearn", "vsearch"]:
                input.append(f"results/qiime2/{short}/queries/{short}/taxonomy_{classifier}.eval.tsv")
    return input

rule eval_all_qiime2:
    input:
        get_all_qiime2_eval

def get_all_qiime2_eval_vsearch(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            for classifier in ["vsearch"]:
                input.append(f"results/qiime2/{short}/queries/{short}/taxonomy_{classifier}.eval.tsv")
    return input

rule eval_all_qiime2_vsearch:
    input:
        get_all_qiime2_eval_vsearch

def get_all_qiime2_eval_sklearn(wildcards):
    input = []
    for key in config["benchmark"].keys():
        for case in cases:
            short = f"{key}.{case.split('-')[0]}"
            for classifier in ["sklearn"]:
                input.append(f"results/qiime2/{short}/queries/{short}/taxonomy_{classifier}.eval.tsv")
    return input

rule eval_all_qiime2_sklearn:
    input:
        get_all_qiime2_eval_sklearn

rule evaluate_qiime2:
    output:
        "results/qiime2/{ref}/queries/{query}/taxonomy_{classifier}.eval.tsv"
    input:
        res="results/qiime2/{ref}/queries/{query}/taxonomy_{classifier}.tsv",
        tax=lambda wildcards: config["qiime2"]["query"][wildcards.query].replace(".fasta", ".tsv"),
    log:
        "results/qiime2/{ref}/queries/{query}/qiime2_{classifier}.eval.log"
    params:
        classifier_str = lambda wildcards: "qiime2_"+wildcards.classifier,
    shell:
        """
        python workflow/scripts/evaluate_classifier.py {input.res} --taxonomy {input.tax} --classifier {params.classifier_str} --output {output} >{log} 2>&1
        """