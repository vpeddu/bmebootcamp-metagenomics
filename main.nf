//Input files 

FASTQ_ch = file("${params.fastq}")

DATABASE_ch = file("${params.krakendb}")

FASTA_ch = Channel
            .fromPath("${params.reference_genomes}**.fasta")
            .map { file -> tuple(file.baseName, file) }



// // adding to kraken database
// yeast_genome_file = file("/Users/vikas/Downloads/GCF_000146045.2_R64_genomic.fna")

// process KrakenDBadd {  
//     //pod secret: "prp-s3-credentials/credentials", mountPath: "/root/.aws/credentials"
//     container 'zdhali/metagenomics:bash'   
//     //pod = [secret: "prp-s3-credentials", mountPath: "/root/.aws"]
//     publishDir "/Users/vikas/Downloads/bootcamp/metagenomics/db_add/"

//      input:
//         file krakendb from DATABASE_ch
//         file yeast_genome from yeast_genome_file

//      output:


//     cpus 6
//     memory 8.Gb 

//     script:
//     """
//     #!/bin/bash

//     echo "adding kraken db"

//     kraken2-build --add-to-library ${yeast_genome} --db ${krakendb}
    
//     echo "downloading taxonomy"

//     kraken2-build --download-taxonomy --db ${krakendb}

//     ls -latr ${krakendb}

//     echo "building"

//     kraken2-build --build --threads ${task.cpus} --db ${krakendb} 
//     """
// }


//Kraken2 run with original database
process Kraken {  

    container 'zdhali/metagenomics:bash'   
    
    publishDir "KrakenOut/"

    input:
        file krakendb from DATABASE_ch
        file fastq from FASTQ_ch

    output:
        tuple file("classified.fastq"), file("unclassified.fastq"), file("kraken_report.tsv"), file("classified.fasta") into classified_ch

    cpus 6
    memory 8.Gb 

    script:
    """
    #!/bin/bash

    echo "kraken running"
    kraken2 --db ${krakendb} --threads ${task.cpus} ${fastq} \
    --classified-out classified.fastq --unclassified-out unclassified.fastq \
    --gzip-compressed --report kraken_report.tsv --use-names

    echo "kraken finished. Logging ls"
    ls -lah

    echo "zipping files"
    #gzip classified.fastq 
    #gzip unclassified.fastq   

    #echo "finding unique taxids"
    #zcat classified.fastq.gz | grep "taxid|"  | cut -d "|" -f2 | sort | uniq > unique_taxids.txt


    cat classified.fastq | awk '{if(NR%4==1) {printf(">%s\\n",substr(\$0,2));} else if(NR%4==2) print;}' > classified.fasta
    """
}


process extract_fastq {
    
    container "quay.io/fhcrc-microbiome/clomp:v0.1.3"
    
    input: 
        tuple file(CLASSIFIED), file(UNCLASSIFIED), file(TSV), file(UNIQUE_TAXIDS) from classified_ch
    output: 
        file('species_taxids.txt')
        file("*.species.fasta") into species_fastq_ch

    script: 
    """
    #!/usr/bin/env python3
import os 
import subprocess
    

newfile = open("species_taxids.txt", "w")
with open('${TSV}') as f:
    for line in f:
        if str(line.split()[3]) == 'S' and int(line.split()[2]) > 50: 
            print(line)
            output = str(line.split()[4])
            name = str(line.split()[5]).replace(" ", "")
            #split_cmd = "zgrep -A4 'taxid|" + output + "' classified.fastq.gz | gzip > " +  name + ".fastq.gz"
            #grep -A 2 -B 1 'AAGTTGATAACGGACTAGCCTTATTTT' file.fq | sed '/--/d'
            #split_cmd = "zgrep -A 2 -B 1 --no-group-separator 'taxid|" + output + "' classified.fastq.gz  | gzip > " +  name + ".fastq.gz"
            split_cmd = "zgrep -A 1 --no-group-separator 'taxid|" + output + "' classified.fasta  > " +  name + ".species.fasta"
            subprocess.call(split_cmd, shell = True)
            print(split_cmd)
            newfile.write(output + "\\n")
newfile.close()



    """


}



// process Assemble {  
//     container 'staphb/flye:2.8'   
    
//     cpus 8 
//     memory 12.Gb

//     input:
//         file (FASTQ) from species_fastq_ch.flatten()
//     script: 

//     """
//     #!/bin/bash


//     echo "logging ls"
//     ls -latr

//     # deleting last line

//     #zcat ${FASTQ} | sed -n '1~4s/^@/>/p;2~4p' > new_fasta.fasta
//     zcat ${FASTQ} | awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' > new_fasta.fasta
//     cat new_fasta.fasta | tr -s ' ' '_' > trimmed.fasta
    
//     number=\$(zgrep ">" new_fasta.fasta | wc -l)
//     echo \$number
//     if [ \$number \\> 100 ]
//     then
//     flye \
//     -t ${task.cpus} \
//     --nano-corr trimmed.fasta \
//     -o filebase.assembled
//     echo "here"
//     fi


    
//     """
// }



process AssembleCanu {  
    container 'shuangbroad/canu:sh_optimize_merylcount_batch'   
    
    cpus 8 
    memory 12.Gb

    input:
        file (FASTA) from species_fastq_ch.flatten()
    script: 

    """
    #!/bin/bash
    

    echo "logging ls"
    ls -latr

    # deleting last line

    #zcat ${FASTA} | sed -n '1~4s/^@/>/p;2~4p' > new_fasta.fasta
    #zcat ${FASTA} | awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' > new_fasta.fasta
    cat ${FASTA} | tr -s ' ' '_' > trimmed.fasta
    
    number=\$(zgrep ">" trimmed.fasta | wc -l)
    echo \$number
    if [ \$number \\> 100 ]
    then
        canu \
        -assemble \
        -p test -d oof \
        genomeSize=5.2m \
        -nanopore-corrected trimmed.fasta
    echo "here"
    fi


    
    """
}
// Running Minimap2 to get alignment
process Minimap2 {  
    pod secret: "prp-s3-credentials/credentials", mountPath: "/root/.aws/credentials"
    container 'nfcore/nanoseq'   
    
    publishDir "Minimap2Out/"

    input:
        tuple val(base), file(reference_fasta) from FASTA_ch
        file fastq from FASTQ_ch

    output: 
        file "${base}.bam"
    cpus 4
    memory 8.Gb 

    script:
    """
    #!/bin/bash

    echo "file ${base}"

    minimap2 -t ${task.cpus} -x map-pb -a ${reference_fasta} ${fastq} | samtools view -Sb -F 4 > ${base}.bam

    """
}
