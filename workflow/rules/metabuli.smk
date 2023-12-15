localrules:
    sintax_fasta_to_gtdbfmt,
    gtdb_to_taxdump,
    refmt_accession2taxid,

rule sintax_fasta_to_gtdbfmt:
    output:
        tax="results/metabuli/{ref}/gtdbfmt/taxonomy.tsv",
        seqs="results/metabuli/{ref}/gtdbfmt/seqs.fasta",
        tbl="results/metabuli/{ref}/gtdbfmt/taxonomy.tbl",
    input:
        lambda wildcards: config["metabuli"]["ref"][wildcards.ref],
    log:
        "results/metabuli/{ref}/gtdbfmt/sintax_fasta_to_gtdbfmt.log"
    params:
        src=srcdir("scripts/coidb_to_gtdbfmt.py"),
        outdir=lambda wildcards, output: os.path.dirname(output.tax),
    shell:
        """
        python {params.src} {input} {params.outdir} > {log} 2>&1
        echo -e "accession" > {output.tbl}
        cut -f1 {output.tax} >> {output.tbl}
        """

rule gtdb_to_taxdump:
    output:
        dump=expand("results/metabuli/{{ref}}/db/taxonomy/{f}.dmp", f = ["delnodes","merged","names","nodes"]),
        tsv="results/metabuli/{ref}/gtdbfmt/taxonomy.taxdump.tsv",
        tbl_wtaxids="results/metabuli/{ref}/gtdbfmt/taxonomy_wTaxIDs.tsv",
        mapfile="results/metabuli/{ref}/gtdbfmt/map.tsv",
    input:
        tax=rules.sintax_fasta_to_gtdbfmt.output.tax,
        tbl=rules.sintax_fasta_to_gtdbfmt.output.tbl,
    log:
        "results/metabuli/{ref}/gtdbfmt/gtdb_to_taxdump.log"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.dump[0]),
    conda:
        "envs/metabuli.yml"
    shell:
        """
        gtdb_to_taxdump.py -t {input.tbl} -o {params.outdir} {input.tax} >{output.tsv} 2>{log}
        cut -f1,2 {output.tbl_wtaxids} > {output.mapfile}
        echo -e "\t|\t\t|" > {params.outdir}/merged.dmp
        """

rule refmt_accession2taxid:
    output:
        a2id="results/metabuli/{ref}/gtdbfmt/accession2taxid.tsv"
    input:
        mapfile=rules.gtdb_to_taxdump.output.mapfile
    run:
        with open(input.mapfile, 'r') as fhin, open(output.a2id, 'w') as fhout:
            for i, line in enumerate(fhin):
                if i ==0:
                    fhout.write("accession\taccession.version\ttaxid\tgi\n")
                    continue
                seqid, taxid = line.rstrip().split("\t")
                fhout.write(f"{seqid[0:-2]}\t{seqid}\t{taxid}\t0\n")       

rule metabuli_add_and_build:
    output:
        expand("results/metabuli/{{ref}}/db/{f}", f=["taxID_list","acc2taxid.map"]),
    input:
        a2id=rules.refmt_accession2taxid.output.a2id,
        seqs=rules.sintax_fasta_to_gtdbfmt.output.seqs,
        dump=expand("results/metabuli/{{ref}}/db/taxonomy/{f}.dmp", f = ["delnodes","merged","names","nodes"]),
    log:
        add="results/metabuli/{ref}/db/add.log",
        build="results/metabuli/{ref}/db/build.log",
    conda: "envs/metabuli.yml"
    params:
        seq_abs = lambda wildcards, input: os.path.abspath(input.seqs),
        dbdir="results/metabuli/{ref}/db",
        local_dbdir="$TMPDIR/metabuli/{ref}/db",
    threads: 20
    resources:
        runtime = 60 * 24 * 5,
        constraint = "mem256GB",
    shell:
        """
        mkdir -p {params.local_dbdir}
        cp -r {params.dbdir}/* {params.local_dbdir}/
        echo {params.seq_abs} > {params.local_dbdir}/list
        metabuli add-to-library {params.local_dbdir}/list {input.a2id} {params.local_dbdir} > {log.add} 2>&1

        find {params.local_dbdir}/library -type f -name '*.fna' > {params.local_dbdir}/library-files.txt
        metabuli build {params.local_dbdir} {params.local_dbdir}/library-files.txt {input.a2id} --threads {threads} >{log.build} 2>&1
        rm -rf {params.local_dbdir}/library
        rm -rf {params.local_dbdir}/taxonomy
        cp {params.local_dbdir}/* {params.dbdir}/
        rm -rf {params.local_dbdir}
        """

rule metabuli_classify:
    output:
        "results/metabuli/{ref}/classify/{run}/{query}_classifications.tsv",
        "results/metabuli/{ref}/classify/{run}/{query}_report.tsv",
    input:
        db=rules.metabuli_add_and_build.output,
        fa=lambda wildcards: config["metabuli"]["query"][wildcards.query],
    log:
        "results/metabuli/{ref}/classify/{run}/{query}.log"
    params:
        dbdir=lambda wildcards, input: os.path.dirname(input.db[0]),
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        jobid=lambda wildcards: wildcards.query,
        options=lambda wildcards: config["metabuli"]["runs"][wildcards.run]["options"]
    conda: "envs/metabuli.yml"
    threads: 20
    resources:
        runtime = 60 * 2,
        constraint="mem256GB"
    shell:
        """
        metabuli classify {input.fa} {params.dbdir} {params.outdir} {params.jobid} --seq-mode 1 --max-ram 256 --threads {threads} {params.options} >{log} 2>&1
        """
