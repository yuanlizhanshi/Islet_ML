import os
import glob


Bowtie2_index = "~/Desktop/hg38/hg38_bowtie2/hg38"

SAMPLES = []

files = glob.glob('rawdata/*.gz')
for i in files:
    SAMPLES.append(os.path.splitext(os.path.splitext(os.path.basename(i))[0])[0])

    
    
rule all:
    input:
        expand("rawdata/{sample}.fastq.gz",sample=SAMPLES),
        expand("clean_fastq/{sample}.fq.gz",sample=SAMPLES),
        expand("rmdup_bam/{sample}_rmdup.bam",sample=SAMPLES),
        expand("rmdup_bam/{sample}_rmdup.bam.bai",sample=SAMPLES),
        expand("bigwig/{sample}_rmdup.bigwig",sample=SAMPLES)


rule QC:
    input:
        "rawdata/{sample}.fastq.gz"
    output:
        "clean_fastq/{sample}.fq.gz"
    threads: 5
    log:
        html = "clean_fastq/{sample}.html",
        json = "clean_fastq/{sample}.json" 
    shell:
        "fastp -w {threads} -i {input} -o {output} --html {log.html} --json {log.json}"

rule Bowtie2_map:
  input:
    "clean_fastq/{sample}.fq.gz"
  output:
    temp("sam/{sample}.sam")
  threads: 10
  log:
    "sam/{sample}_mapping_log.txt"
  shell:
    "bowtie2 -p {threads} -q {input} -x {Bowtie2_index} -S {output} 2>{log}"


rule samtools_view:
  input:
    "sam/{sample}.sam"
  output:
    temp('sam/{sample}.bam')
  threads: 10
  shell:
    'samtools view -q 10 -bh -o {output} {input}'


rule samtools_sort:
  input:
    'sam/{sample}.bam'
  output:
    temp('sortedbam/{sample}.bam')
  threads: 10
  shell:
    'samtools sort -@ {threads} -o {output} {input}'

rule samtools_remove_duplication:
  input:
    'sortedbam/{sample}.bam'
  output:
    bam = 'rmdup_bam/{sample}_rmdup.bam',
  shell:
    'samtools rmdup {input} {output.bam}'

rule bam_index:
  input:
    'rmdup_bam/{sample}_rmdup.bam'
  output:
    'rmdup_bam/{sample}_rmdup.bam.bai'
  shell:
    "samtools index {input}"


rule bam_to_bigwig:
  input:
    bam = 'rmdup_bam/{sample}_rmdup.bam',
    bam_index = 'rmdup_bam/{sample}_rmdup.bam.bai'
  output:
    'bigwig/{sample}_rmdup.bigwig'
  threads: 10
  log:
    'bigwig/{sample}.bigwig.log'
  shell:
    'bamCoverage -p {threads} -bs 20 -b {input.bam} -o {output} --normalizeUsing CPM 2>{log}'

