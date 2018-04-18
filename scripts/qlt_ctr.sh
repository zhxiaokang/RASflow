# use FastQC to check the quality of reads
KEY=/export/jonassenfs/xiaokangz/project/data/RNA-Seq-analysis/key/u11bd@nelstor0.cbu.uib.no.txt
FASTQC=/Home/ii/xiaokangz/tools/bioinfo/FastQC/fastqc
FILEID=/export/jonassenfs/xiaokangz/project/RNA-Seq-analysis/configs/fastq.file.names.id
READSDIR=/export/jonassenfs/xiaokangz/project/data/RNA-Seq-analysis/reads
OUTDIR=/export/jonassenfs/xiaokangz/project/RNA-Seq-analysis/output/fastqc
NELSIN=u11bd@nelstor0.cbu.uib.no:Projects/UiB_Goksoyr_dCod_1_2017/dCod_Batch3
NELSOUT=u11bd@nelstor0.cbu.uib.no:Personal/UiB_Goksoyr_dCod_1_2017/dCod_Batch3/Essa/fastqc
NCORE=10

while IFS='' read -r i || [[ -n "$i" ]]; do
    echo Start processing sample $i at: `date`
    echo Copy paired reads from NeLS to server......
    forward=${i}_S*_R1_001.fastq.gz
    reverse=${i}_S*_R2_001.fastq.gz
    scp -i $KEY ${NELSIN}/$forward $READSDIR
    scp -i $KEY ${NELSIN}/$reverse $READSDIR
	
	##  ===================  Use FastQC to check the quality of reads  =====================
	
	echo Start analysing the quality at: `date`
	$FASTQC -o $OUTDIR $READSDIR/$forward $READSDIR/$reverse
	echo Finish analysing the quality at: `date`

	##====================  Deal with the files  =====================

	# remove the reads files
	rm $READSDIR/$forward $READSDIR/$reverse

	echo Finish processing sample $i at: `date`
done < $FILEID

# Send local output files to NeLS
echo Sending files to NeLS......
    # send the quality control report to NeLS
scp -i $KEY $OUTDIR/* $NELSOUT


