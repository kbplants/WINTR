export TMPDIR=/mnt/winterprojectceph/nf_tmp
nextflow run main_qc.nf -resume -w /mnt/winterprojectceph/work -with-report report401_600.html -with-trace trace401_600.txt > log401_600.txt
