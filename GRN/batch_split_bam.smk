import glob
import os

directory_path = "raw_bam"


bam_files = []
files = glob.glob(directory_path + '/*.bam')
for i in files:
    bam_files.append(os.path.splitext(os.path.basename(i))[0])

rule all:
  input:
    expand("raw_bam/{sample}.bam",sample=bam_files),
    expand("raw_bam/{sample}.bam.bai",sample=bam_files),
    expand("split_bam/{sample}.report.txt",sample=bam_files)

rule index:
    input:
        "raw_bam/{sample}.bam"
    output:
        'raw_bam/{sample}.bam.bai'
    shell:
        "samtools index {input}"
        
rule split_bam:
    input:
        bam = "raw_bam/{sample}.bam",
        bai = "raw_bam/{sample}.bam.bai",
        meta = "metadata/{sample}.tsv"
    output:
        'split_bam/{sample}.report.txt'
    shell:
        "python SplitBamCellTypes.py --bam {input.bam} --meta {input.meta} --id {wildcards.sample} --n_trim 5 --min_MQ 30  --outdir split_bam"        
        
        
        