SAMPLES, = glob_wildcards("/data/{sample}.gbk")
   
rule all:
	input: "/data/result/README.md"
# expand("/data/classification/{sample}_assigned_locus.tsv", sample=SAMPLES)
   
rule aa_seq:
    input: "/data/{sample}.gbk"
    output: temp("/data/aa_fasta/{sample}_aa.fasta")
    shell: "python -W ignore /STEC_prediction/script/gbk_aaseq_extraction.py -i {input}"    
     
rule locus:
    input: "/data/{sample}.gbk"
    output:
        temp("/data/locus/{sample}.tsv")
    shell: "python /STEC_prediction/script/locus_tag.py -i {input}"          
 
rule diamond:
    input: "/data/aa_fasta/{sample}_aa.fasta"
    output: 
        blastp = temp("/data/blastp/blastp/{sample}.blastp"),
        blastp_stx = temp("/data/blastp/blastp_stx/{sample}_stx.blastp")
    run:
        shell("mkdir -p /data/blastp && diamond blastp -q {input} -d /STEC_prediction/reference.dmnd -o {output.blastp} -f 6")
        shell("diamond blastp -q {input} -d /STEC_prediction/stx.dmnd -o {output.blastp_stx} -f 6")
    
rule pirate:
    input: 
        diamond="/data/blastp/blastp/{sample}.blastp",
        locus="/data/locus/{sample}.tsv" 
    output:temp("/data/classification/{sample}_assigned_locus.tsv")
    shell: "python /STEC_prediction/script/pirate_gene_classification_locus.py -d {input.diamond} -l {input.locus}"

rule prediction:
    input: 
        expand("/data/classification/{sample}_assigned_locus.tsv", sample=SAMPLES), 
        expand("/data/blastp/blastp/{sample}.blastp", sample=SAMPLES),
        expand("/data/blastp/blastp_stx/{sample}_stx.blastp", sample=SAMPLES)
    output: "/data/result/README.md"
    run:
        shell("python /STEC_prediction/script/prediction.py")
        shell("echo 'If you want to add other .gbk files and run it again, please move the data in [result] folder and remove it. Then run the pipeline again.' > /data/result/README.md")
