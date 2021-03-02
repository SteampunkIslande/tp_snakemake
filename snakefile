# bwa index genom/ecoli.fa 
# bwa mem genom/ecoli.fa sample1.fastq > sample1.sam
# samtools sort -O BAM sample1.sam > sample1.bam

#bwa mem -R "RG\tID:sample1\tSM:sample1"

import glob
SAMPLES = [ f.split('.')[0] for f in glob.glob('*.fastq')]

rule align:
   input:
      "genom/ecoli.fa",
      "{sample}.fastq"
   output:
      "{sample}.sam"
   shell:
      "bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' {input} > {output}"

rule sam2bam:
   input:
      "{sample}.sam"
   output:
      "{sample}.bam"
   shell:
      "samtools sort -O BAM {input} > {output} ; samtools index {output}"

rule bam2vcf:
   input:
      "genom/ecoli.fa",
      "{sample}.bam"
   output:
      "{sample}.vcf"
   shell:
      "freebayes -f {input} > {output}"

rule vcf2bgz:
   input:
      "{sample}.vcf"
   output:
      "{sample}.vcf.gz"
   shell:
      "bgzip {input} ; tabix {output}"

rule mergeAll : 
    input:
        expand("{sample}.vcf.gz", sample=SAMPLES)
    output:
        "all.vcf.gz"
    shell:
        "bcftools merge {input}|bgzip> {output}"
