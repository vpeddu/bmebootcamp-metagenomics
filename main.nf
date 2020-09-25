
// process Testing {  
//     //pod secret: "prp-s3-credentials/credentials", mountPath: "/root/.aws/credentials"
//     container 'tbattaglia/kraken2'   
//     pod = [secret: "prp-s3-credentials", mountPath: "/root/.aws"]
//     publishDir "s3://bmebootcamp/vpeddu/ahhhhhhhhhhhhhh/"

//     // input:
//     //     file krakendb from DATABASE_ch
//     //     file fastq from FASTQ_ch

//      output:
//          file "credentials.txt"

//     cpus 2
//     memory 7.Gb 

//     script:
//     """
//     #!/bin/bash

//     cat /root/.aws/credentials > credentials.txt


//     """
// }

 //FASTQ_ch = file("${params.fastq}")
 DATABASE_ch = file("${params.krakendb}")

FASTA_ch = Channel
            .fromPath("${params.reference_genomes}**.fasta")
            .map { file -> tuple(file.baseName, file) }


// process Kraken {  
//     pod secret: "prp-s3-credentials/credentials", mountPath: "/root/.aws/credentials"
//     container 'tbattaglia/kraken2'   
    
//     publishDir "KrakenOut/"

//     input:
//         file krakendb from DATABASE_ch
//         file fastq from FASTQ_ch

//     output:
//         file "classified.fastq"
//         file "unclassified.fastq"
//         file "kraken_report.tsv"

//     cpus 2
//     memory 7.Gb 

//     script:
//     """
//     #!/bin/bash

//     echo "kraken running"
//     kraken2 --db ${krakendb} --threads ${task.cpus} ${fastq} --quick \
//     --classified-out classified.fastq --unclassified-out unclassified.fastq \
//     --gzip-compressed --report kraken_report.tsv --use-names

//     echo "kraken finished. Logging ls"
//     ls -latr


//     """
// }

// process Minimap2 {  
//     pod secret: "prp-s3-credentials/credentials", mountPath: "/root/.aws/credentials"
//     container 'biocontainers/minimap2:v2.15dfsg-1-deb_cv1'   
    
//     publishDir "Minimap2Out/"

//     input:
//         tuple val(base), file(reference_fasta) from FASTA_ch
//         file fastq from FASTQ_ch

//     output: 
//         file "${base}.sam"
//     cpus 2
//     memory 4.Gb 

//     script:
//     """
//     #!/bin/bash

//     echo "file ${base}"

//     minimap2 -t ${task.cpus} -x map-pb -a ${reference_fasta} ${fastq} > ${base}.sam

//     """
// }


yeast_genome_file = file("/Users/vikas/Downloads/GCF_000146045.2_R64_genomic.fna")

process KrakenDBadd {  
    //pod secret: "prp-s3-credentials/credentials", mountPath: "/root/.aws/credentials"
    container 'tbattaglia/kraken2'   
    //pod = [secret: "prp-s3-credentials", mountPath: "/root/.aws"]
    publishDir "/Users/vikas/Downloads/bootcamp/metagenomics/db_add/"

     input:
        file krakendb from DATABASE_ch
        file yeast_genome from yeast_genome_file

     output:


    cpus 4
    memory 7.Gb 

    script:
    """
    #!/bin/bash

    # remove one zafrins dockerfile finishes building
    apt install -y rsync

    echo "adding kraken db"

    kraken2-build --add-to-library ${yeast_genome} --db ${krakendb}
    
    echo "download taxonomy"

    kraken2-build --download-taxonomy --db ${krakendb}

    ls -latr ${krakendb}

    echo "building"

    kraken2-build --build --threads ${task.cpus} --db ${krakendb} 

    """
}