def mem_allowed(wildcards, threads):
    return max(threads * 6400, 6400)

rule all:
    input:
        expand("results/alignment/{db}/clustalo-aligned.fasta", db=config["database"].keys()),

rule clustalomega_create_profile:
    """
    Create a profile alignment of the database sequences
    """
    output:
        aln="results/alignment/{db}/clustalo-profile.fasta",
        filelist="results/alignment/{db}/fileList.txt"
    input:
        fasta = lambda wildcards: config["database"][wildcards.db]["fasta"],
        taxfile = lambda wildcards: config["database"][wildcards.db]["taxonomy"],
    log:
        "results/alignment/{db}/clustalomega_create_profile.log"
    threads: 20
    envmodules:
        "bioinfo-tools",
        "clustalo/1.2.4",
        "biopython/1.76"
    conda:
        "../envs/clustalo.yml"
    resources:
        runtime = 60 * 24 * 2, 
        mem_mb = mem_allowed,
    params:
        rank="family",
        outdir=lambda wildcards, output: os.path.dirname(output.aln)
    shell:
        """
        python workflow/scripts/clustalo-rank-profile.py \
            --taxfile {input.taxfile} --seqsfile {input.fasta} \
            --profile {output.aln} --outdir {params.outdir} --rank {params.rank} --threads {threads} > {log} 2>&1
        """

rule clustalomega_align_to_profile:
    output:
        "results/alignment/{db}/_aligned/file_{i}.aln"
    input:
        profile=rules.clustalomega_create_profile.output.aln,
        fasta="results/alignment/{db}/_files/file_{i}.fna"
    log:
        "benchmark/{db}/_logs/file_{i}.log"
    envmodules:
        "bioinfo-tools",
        "clustalo/1.2.4",
    conda:
        "../envs/clustalo.yml"
    threads: 20
    resources:
        runtime = 60 * 12, 
        mem_mb = mem_allowed,
    shell:
        """
        clustalo -i {input.fasta} -o {output} --outfmt=a2m -t DNA --profile1 {input.profile} --threads {threads} > {log} 2>&1
        """

def get_files(wildcards):
    """
    Return a list of files to concatenate
    """
    files = []
    try:
        if os.path.exists(config["alignment"]["fileList"]):
            with open(config["fileList"], 'r') as fhin:
                for line in fhin:
                    f = line.strip().replace(".fna", ".aln")
                    files.append(f)
        return files
    except KeyError:
        return []

rule clustalomega_concat_alns:
    output:
        "results/alignment/{db}/clustalo-aligned.fasta"
    input:
        get_files
    log:
        "results/alignment/{db}/clustalo-concat.log"
    envmodules:
        "bioinfo-tools",
        "clustalo/1.2.4",
    conda:
        "../envs/clustalo.yml"
    threads: 20
    resources:
        runtime = 60 * 12, 
        mem_mb = mem_allowed,
    shell:
        """
        python workflow/scripts/concat-aln.py --fasta {input} -o {output} > {log} 2>&1
        """