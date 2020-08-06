####
# InSilicoSCR for sc2 guides pipeline
# author: Chunyu Zhao chunyu.zhao@czbiohub.org
# time: 2020-08-05
####

import os
import sys
import configparser
from Bio import SeqIO


GUIDES = [record.id for record in SeqIO.parse(config["guide_fasta"], "fasta")]

rule _all_neighbors:
    input:
        expand(config["project_dir"] + "/1_neighbors/{guide}_hd.4.txt", guide=GUIDES)

rule generate_neighbors:
    input:
        config["guide_fasta"]
    output:
        config["project_dir"] + "/1_neighbors/{guide}_hd.4.txt"
    shell:
        "python isscrlib/gen_neighbors.py {input} {output}"


rule _all_scr_db:
    input:
        expand(config["project_dir"] + "/2_kq_vn/{guide}/classified_reads.tsv", guide=GUIDES),
        expand(config["project_dir"] + "/2_kq_rb/{guide}/classified_reads.tsv", guide=GUIDES)

rule scr_viral_neighbors:
    input:
        db = config["krakenuniq"]["viral_neighbors_db"],
        neighbor = config["project_dir"] + "/1_neighbors/{guide}_hd.4.txt"
    output:
        config["project_dir"] + "/2_kq_vn/{guide}/classified_reads.tsv"
    threads: 8
    params:
        config["project_dir"] + "/2_kq_vn/{guide}"
    shell:
        """
        krakenuniq --fasta-input {input.neighbor}  \
            --classified-out {params}/classified_out.fasta \
            --db {input.db} --threads {threads} \
            --report-file {params}/report_file.tsv \
            --output {params}/read_classification.tsv

        grep C {params}/read_classification.tsv > {params}/classified_reads.tsv
        """

rule scr_refseq_bacteria:
    input:
        db = config["krakenuniq"]["refseq_bacteir_db"],
        neighbor = config["project_dir"] + "/1_neighbors/{guide}_hd.4.txt"
    output:
        config["project_dir"] + "/2_kq_rb/{guide}/classified_reads.tsv"
    threads: 8
    params:
        "2_kq_rb/{guide}"
    shell:
        """
        krakenuniq --fasta-input {input.neighbor}  \
            --classified-out {params}/classified_out.fasta \
            --db {input.db} --threads {threads} \
            --report-file {params}/report_file.tsv \
            --output {params}/read_classification.tsv

        grep C {params}/read_classification.tsv > {params}/classified_reads.tsv
        """

rule _all_db_classified_kmers:
    input:
        expand(config["project_dir"] + "/3_classified_neighbors_vn/{guide}.tsv", guide=GUIDES),
        expand(config["project_dir"] + "/3_classified_neighbors_rb/{guide}.tsv", guide=GUIDES)

rule classify_neighbors_vn:
    input:
        config["project_dir"] + "/2_kq_vn/{guide}/classified_reads.tsv"
    output:
        config["project_dir"] + "/3_classified_neighbors_vn/{guide}.tsv"
    params:
        guides = config["guide_fasta"],
        db = config["krakenuniq"]["viral_neighbors_db"],
    shell:
        "Rscript /mnt/chunyu_6TB/cas13/sc2-guide-InSilicoSCR/isscrlib/classify_neighbors.R {wildcards.guide} {input} {output} {params.db} {params.guides}"

rule classify_neighbors_rb:
    input:
        config["project_dir"] + "/2_kq_rb/{guide}/classified_reads.tsv"
    output:
        config["project_dir"] + "/3_classified_neighbors_rb/{guide}.tsv"
    params:
        guides = config["guide_fasta"],
        db = config["krakenuniq"]["viral_neighbors_db"],
    shell:
        "Rscript isscrlib/classify_neighbors.R {wildcards.guide} {input} {output} {params.db} {params.guides}"
