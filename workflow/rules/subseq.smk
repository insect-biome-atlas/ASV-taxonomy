import os

localrules: 
    gather_subseqs,
    translate_db,
    download_hmm,
    extract_hmm_subseq,
    trim_subseq_aln,
    gather_subseqs


splits=[f'{x:03d}' for x in list(range(1,10001))]

rule subseq:
    input:
        expand("benchmark/{db}/_trimmed/{split}.fasta", split=splits, db=config["benchmark"].keys()),


rule translate_db:
    """
    Translate each reference sequence in all six reading frames and choose the longest
    """
    output:
        faa = "benchmark/{db}/db.faa",
        coord = "benchmark/{db}/db.coord"
    input:
        fasta = lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
    log:
        "benchmark/{db}/translate_db.log"
    shell:
        """
        python workflow/scripts/translate.py --coord {output.coord} {input.fasta} >{output.faa} 2>{log}
        """

rule split_translations:
    """
    Split translated aa sequence file tino 10,000 chunks
    """
    output:
        expand("benchmark/{{db}}/_splits/split{split}.faa", split=splits)
    input:
        rules.translate_db.output.faa,
    log:
        "benchmark/{db}/split_translations.log"
    conda:
        "../envs/seqkit.yml"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    resources:
        runtime = 120,
    threads: 4
    shell:
        """
        seqkit split2 -O {params.outdir} -j {threads} -p 10000 --by-part-prefix split -1 {input} >{log} 2>&1
        """

rule download_hmm:
    """
    Download COX1 HMM from Interpro
    """
    output:
        "data/PF00115.hmm"
    shell:
        """
        curl -L  "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00115?annotation=hmm" | gunzip -c > {output}
        """

rule hmmsearch:
    output:
        "benchmark/{db}/_hmmsearch/{split}.out"
    input:
        split="benchmark/{db}/_splits/split{split}.faa",
        hmm="data/PF00115.hmm"
    log:
        "benchmark/{db}/_hmmsearch/{split}.log"
    conda:
        "../envs/hmmer.yml"
    threads: 4
    resources:
        runtime = 60,
    shell:
        """
        hmmsearch --domtblout {output} --cut_tc --cpu {threads} {input.hmm} {input.split} > /dev/null 2>{log}
        """

rule extract_hmm_subseq:
    output:
        "benchmark/{db}/_subseq/{split}.fasta"
    input:
        hmmout="benchmark/{db}/_hmmsearch/{split}.out",
        dnafile=lambda wildcards: config["benchmark"][wildcards.db]["fasta"],
        coord="benchmark/{db}/db.coord"
    log:
        "benchmark/{db}/_subseq/{split}.log"
    params:
        min_len=lambda wildcards: config["benchmark"][wildcards.db]["subseq_min_len"],
        hmm_from=lambda wildcards: config["benchmark"][wildcards.db]["hmm_from"],
        hmm_to=lambda wildcards: config["benchmark"][wildcards.db]["hmm_to"],
    shell:
        """
        python workflow/scripts/extract_hmm_subseq.py -m {params.min_len} --hmm_from {params.hmm_from} \
            --hmm_to {params.hmm_to} {input.hmmout} {input.coord} {input.dnafile} > {output} 2>{log}
        """

def get_asv_seq(wildcards):
    if "asv_seq" in config["benchmark"][wildcards.db].keys():
        return ">ASV\n" + config["benchmark"][wildcards.db]["asv_seq"] + "\n"
    else:
        return ""

rule align_subseqs:
    output:
        "benchmark/{db}/_aln/{split}.aln"
    input:
        "benchmark/{db}/_subseq/{split}.fasta"
    log:
        "benchmark/{db}/_aln/{split}.log"
    conda:
        "../envs/clustalo.yml"
    threads: 4
    resources:
        runtime = 30,
    params:
        asv_seq = lambda wildcards: get_asv_seq(wildcards),
        tmpdir = "$TMPDIR/clustalo.{db}",
        out = "$TMPDIR/clustalo.{db}/subseqs.aln"
    shell:
        """
        mkdir -p {params.tmpdir}
        echo -e "{params.asv_seq}" > {params.tmpdir}/input.fasta
        cat {input} >> {params.tmpdir}/input.fasta
        clustalo -i {params.tmpdir}/input.fasta -o {params.out} --outfmt=a2m -t DNA --threads={threads} -v 2>{log}
        mv {params.out} {output}
        rm -rf {params.tmpdir}
        """

rule trim_subseq_aln:
    """
    Trim alignment to match only the columns of the ASV sequence in the alignment
    """
    output:
        "benchmark/{db}/_trimmed/{split}.fasta"
    input:
        "benchmark/{db}/_aln/{split}.aln"
    log:
        "benchmark/{db}/_aln/{split}.log"
    params:
        min_len=300,
    shell:
        """
        python workflow/scripts/trimaln.py -m {params.min_len} {input} > {output} 2>{log}
        """

rule gather_subseqs:
    output:
        "benchmark/{db}/subseqs.fasta"
    input:
        expand("benchmark/{{db}}/_trimmed/{split}.fasta", split=splits),
    log:
        "benchmark/{db}/gather_subseqs.log"
    shell:
        """
        cat {input} > {output} 2>{log}
        """

