nextflow.enable.dsl=21


ctab_ch = Channel.fromPath("${params.data_dir}/stringtie_out/*/*_t_data_with_tpm.ctab")
workflow {

    fastq_ch = Channel.fromPath("${params.data_dir}/*.fastq.gz")
        .map { fastq ->
            def sample = fastq.getBaseName().replaceFirst(/\.filter-RNA.fastq$/, '')
            tuple(sample, fastq)
        }

    aligned_bam_ch = fastq_ch | split_reads | fastqc_before | trim_reads | fastqc_after | align_reads
    gtf_ctab_ch = aligned_bam_ch | run_stringtie
    
    gtf_ch = gtf_ctab_ch.map { sample, gtf, _ -> gtf }
    ctab_ch = gtf_ctab_ch.map { sample, _, ctab -> ctab }


    gtf_ch.collect() | generate_sample_list | run_prepde
    ctab_ch.collect() | generate_ctab_sample_list | build_tpm_matrices
}



process split_reads {
    input:
    tuple val(sample), path(fastq)
    output:
    tuple val(sample), path("${sample}_R1.fastq"), path("${sample}_R2.fastq")

    script:
    """
    echo "Splitting reads for sample: ${sample}"
    reformat.sh in=${fastq} out1=${sample}_R1.fastq out2=${sample}_R2.fastq
    """
}

process fastqc_before {
    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path(R1), path(R2)

    script:
    """
    mkdir -p ${params.data_dir}/fastqc_before
    fastqc ${R1} ${R2} -o ${params.data_dir}/fastqc_before
    """
}

process trim_reads {
    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path("${sample}_R1_paired.fastq"), path("${sample}_R2_paired.fastq")

    script:
    """
    trimmomatic PE -threads ${params.threads} \
        ${R1} ${R2} \
        ${sample}_R1_paired.fastq ${sample}_R1_unpaired.fastq \
        ${sample}_R2_paired.fastq ${sample}_R2_unpaired.fastq \
        ILLUMINACLIP:${params.trim_db}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    mkdir -p ${params.data_dir}/trimmed
    cp ${sample}_R1_paired.fastq ${params.data_dir}/trimmed/
    cp ${sample}_R2_paired.fastq ${params.data_dir}/trimmed/
    """
}

process fastqc_after {
    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path(R1), path(R2)

    script:
    """
    mkdir -p ${params.data_dir}/fastqc_after_trim
    fastqc ${R1} ${R2} -o ${params.data_dir}/fastqc_after_trim
    """
}

process align_reads {
    input:
    tuple val(sample), path(R1), path(R2)
    output:
    tuple val(sample), path("${sample}.sorted.bam")

    script:
    """
    sam=${sample}.sam
    bam=${sample}.sorted.bam

    hisat2 -p ${params.threads} -x ${params.hisat_index}/genome -1 ${R1} -2 ${R2} -S \$sam
    samtools view -bS \$sam | samtools sort -o \$bam
    samtools index \$bam
    rm \$sam
    """
}



process run_stringtie {
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



