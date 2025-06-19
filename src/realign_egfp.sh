egfp_fasta=$1
reads_folder=$2

echo `date` "processing pool A"
bbmap.sh ref=${egfp_fasta} in=${reads_folder}/P3531_SP255_003_PoolA_S3_L001_R2_001.fastq.gz threads=32 mappedonly=t scafstats=statsA.txt out=outA.sam k=8 minid=0.2 local=t basecov=covA.txt
echo `date` "processing pool B"
bbmap.sh ref=${egfp_fasta} in=${reads_folder}/P3531_SP255_006_PoolB_S4_L001_R2_001.fastq.gz threads=32 mappedonly=t scafstats=statsB.txt out=outB.sam k=8 minid=0.2 local=t basecov=covB.txt
echo `date` "processing pool C"
bbmap.sh ref=${egfp_fasta} in=${reads_folder}/P3531_SP255_009_PoolC_S5_L001_R2_001.fastq.gz threads=32 mappedonly=t scafstats=statsC.txt out=outC.sam k=8 minid=0.2 local=t basecov=covC.txt
