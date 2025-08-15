nextflow.enable.dsl=2


workflow {
        samples_ch = Channel.fromPath('sample401_600.csv')
    .flatMap { file -> file.readLines() }
    .map { line ->
        def filename = line.trim()

        // Define both possible paths
        def path1 = file("/mnt/winterprojectceph/wintr_jgi_samples/data/${filename}")
        def path2 = file("/mnt/winterprojectceph/wintr_fastq/${filename}")

        // Choose the first existing file, else null
        def fastq = path1.exists() ? path1 : (path2.exists() ? path2 : null)

        if (fastq) {
            def sample = fastq.getBaseName().replaceFirst(/\.filter-RNA\.fastq\.gz$/, '')
            tuple(sample, fastq)
        } else {
            // Returning null here will filter out later
            null
        }
    }
    .filter { it != null }


    split_ch = samples_ch | split_reads
    fastqc_before_ch = split_ch | fastqc_before

    fastqc_before_ch.collect() | multiqc_before

    trimmed_ch = fastqc_before_ch | trim_reads | fastqc_after

    trimmed_ch.collect() | multiqc_after
    
    aligned_bam_ch = trimmed_ch | align_reads
    gtf_ctab_ch = aligned_bam_ch | run_stringtie
    
    gtf_ch = gtf_ctab_ch.map { sample, gtf, _ -> gtf }
    ctab_ch = gtf_ctab_ch.map { sample, _, ctab -> ctab }


    gtf_ch.collect() | generate_sample_list | run_prepde
    ctab_ch.collect() | generate_ctab_sample_list | build_tpm_matrices
}




process split_reads {
    maxForks 20
    scratch true
    input:
    tuple val(sample), path(fastq)
    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz")

    errorStrategy 'ignore'

    script:
    """
    echo "Checking integrity of ${fastq}..."
    if gzip -t ${fastq} 2>/dev/null; then
        echo "Splitting reads for sample: ${sample}"
        reformat.sh in=${fastq} out1=${sample}_R1.fastq.gz out2=${sample}_R2.fastq.gz tossbrokenreads=t

        if [ -s ${sample}_R1.fastq.gz ] && [ -s ${sample}_R2.fastq.gz ]; then
            echo "Split successful."
        else
            echo "WARNING: Output files are empty. Skipping sample ${sample}."
            echo "${sample}" >> skipped_samples.txt
            exit 1
        fi
    else
        echo "WARNING: Corrupted FASTQ file detected for sample ${sample}, skipping..."
        echo "${sample}" >> skipped_samples.txt
        exit 1
    fi

    """
}

process fastqc_before {
    maxForks 20
    scratch true
    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path(R1), path(R2)
    
    publishDir "${params.data_dir}/fastqc_before", mode: 'copy'

    script:
    """
    mkdir -p ${params.data_dir}/fastqc_before
    fastqc ${R1} ${R2} -o ${params.data_dir}/fastqc_before
    """
}

process multiqc_before {
    maxForks 20
    scratch true
    input:
    val sample

    output:
    path "multiqc_before"
    
    publishDir "${params.data_dir}/multiqc_before", mode: 'copy'

    script:
    """
    mkdir -p multiqc_before
    multiqc ${params.data_dir}/fastqc_before -o multiqc_before

    """
}


process trim_reads {
    maxForks 20
    scratch true

    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path("${sample}_R1_paired.fastq.gz"), path("${sample}_R2_paired.fastq.gz")
   

    script:
    """
    trimmomatic PE -threads ${params.threads} \
        ${R1} ${R2} \
        ${sample}_R1_paired.fastq.gz ${sample}_R1_unpaired.fastq.gz \
        ${sample}_R2_paired.fastq.gz ${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${params.trim_db}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        -phred33

    #mkdir -p ${params.data_dir}/trimmed
    #cp ${sample}_R1_paired.fastq.gz ${params.data_dir}/trimmed/
    #cp ${sample}_R2_paired.fastq.gz ${params.data_dir}/trimmed/
    """
}


process fastqc_after {
    maxForks 20
    scratch true

    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path(R1), path(R2)

    publishDir "${params.data_dir}/fastqc_after_trim", mode: 'copy'

    script:
    """
    mkdir -p ${params.data_dir}/fastqc_after_trim
    fastqc ${R1} ${R2} -o ${params.data_dir}/fastqc_after_trim
    """
}

process multiqc_after {
    maxForks 20
    scratch true

    input:
    val sample

    output:
    path "multiqc_after"
    
    publishDir "${params.data_dir}/multiqc_after", mode: 'copy'

    script:
    """
    mkdir -p multiqc_after
    multiqc ${params.data_dir}/fastqc_after_trim -o multiqc_after
    """
}


process align_reads {
    maxForks 15
    scratch true 
  
    input:
    tuple val(sample), path(R1), path(R2)

    output:
    tuple val(sample), path("${sample}.sorted.bam")

    script:
    """
    STAR \\
        --runThreadN ${params.threads} \\
        --genomeDir ${params.star_index} \\
        --readFilesIn ${R1} ${R2} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample}_ \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM 25000000000

    mv ${sample}_Aligned.sortedByCoord.out.bam ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
   
    mkdir -p ${params.data_dir}/star_logs
    
    cp ${sample}_Log.out ${params.data_dir}/star_logs/
    cp ${sample}_Log.final.out ${params.data_dir}/star_logs/
    cp ${sample}_Log.progress.out ${params.data_dir}/star_logs/

    samtools flagstat ${sample}.sorted.bam > ${params.data_dir}/star_logs/${sample}_flagstat.txt

    """
}


process run_stringtie {
    maxForks 20
    scratch true
    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.gtf"), path("${sample}_t_data_with_tpm.ctab")

    script:
    """
    mkdir -p ${params.data_dir}/stringtie_out/${sample}
    stringtie ${bam} \\
        -G ${params.gtf_file} \\
        -e -B -p ${params.threads} \\
        -o ${params.data_dir}/stringtie_out/${sample}/${sample}.gtf

    python3 ${params.add_tpm_script} \\
        --gtf ${params.data_dir}/stringtie_out/${sample}/${sample}.gtf \\
        --ctab ${params.data_dir}/stringtie_out/${sample}/t_data.ctab \\
        --output ${params.data_dir}/stringtie_out/${sample}/${sample}_t_data_with_tpm.ctab

    cp ${params.data_dir}/stringtie_out/${sample}/${sample}.gtf ${sample}.gtf
    cp ${params.data_dir}/stringtie_out/${sample}/${sample}_t_data_with_tpm.ctab ${sample}_t_data_with_tpm.ctab
    """
}



process generate_ctab_sample_list {
    maxForks 20
    scratch true
    input:
    path(ctab_files)

    output:
    path("ctab_sample_lst.txt")

    script:
    """
    mkdir -p ${params.data_dir}/result_df
    > ctab_sample_lst.txt
    for ctab in ${ctab_files}; do
        sample=\$(basename \$ctab _t_data_with_tpm.ctab)
        echo -e "\$sample\\t\$(realpath \$ctab)" >> ctab_sample_lst.txt
    done

    cp ctab_sample_lst.txt ${params.data_dir}/result_df/
    """
}



process generate_sample_list {
    maxForks 20
    scratch true
    input:
    path(gtfs)

    output:
    path("sample_lst.txt")

    script:
    """
    mkdir -p ${params.data_dir}/result_df
    > sample_lst.txt
    for gtf in ${gtfs}; do
        sample=\$(basename \$gtf .gtf)
        echo -e "\$sample\\t\$(realpath \$gtf)" >> sample_lst.txt
    done
    cp sample_lst.txt ${params.data_dir}/result_df/
    """
}

process run_prepde {
    maxForks 20
    scratch true
    input:
    path(sample_list)

    output:
    path("gene_count_matrix.csv")
    path("transcript_count_matrix.csv")

    script:
    """
    if [ ! -s ${sample_list} ]; then
        echo " ERROR: sample_lst.txt is empty!" >&2
        exit 1
    fi

    echo "Contents of sample list:"
    cat ${sample_list}

    python3 ${params.prepde_script} -i ${sample_list}
    cp gene_count_matrix.csv transcript_count_matrix.csv ${params.data_dir}/result_df/
    """
}



process build_tpm_matrices {
    maxForks 20
    scratch true
    input:
    path(ctab_sample_list)

    output:
    path("tpm_transcript_matrix.tsv")
    path("tpm_gene_matrix.tsv")

    script:
    """
    if [ ! -s ${ctab_sample_list} ]; then
        echo "ERROR: ctab_sample_lst.txt is empty!" >&2
        exit 1
    fi

    echo "Contents of ctab sample list:"
    cat ${ctab_sample_list}

    python3 ${params.build_tpm_script} \\
      --sample_list ${ctab_sample_list} \\
      --out_prefix tpm

    cp tpm_transcript_matrix.tsv tpm_gene_matrix.tsv ${params.data_dir}/result_df/
    """
}


