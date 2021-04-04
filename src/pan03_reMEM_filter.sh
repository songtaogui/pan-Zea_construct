#!/usr/bin/env bash
if [ $# -ne 6 ]; then echo "usage $0 <pre> <unrefseq> <MemRef> <coverage> <identity> <cpu> " >&2; exit 1;fi

pre=$1
unalnseq=$2
bwaidx=$3
export tcvrg=$4
export tpid=$5
cpu=$6

if [ -s ${pre}.unaligned.MEMfiltered.fa.gz ];then
	echo "Skip running. Output exists: ${pre}.unaligned.MEMfiltered.fa.gz" >&2
	seqkit stats -j $cpu $unalnseq ${pre}.unaligned.MEMfiltered.fa.gz >&2
	exit
fi

if [ -s ${pre}_reMEM.paf.gz ];then
	echo "Skip running MEM. Output file exist: ${pre}_reMEM.paf.gz " >&2
else
	echo "Aligning unaligned seq back to REF with bwa MEM ..." >&2
	# reMEM to REF and filter
	bwa mem -t $cpu $bwaidx $unalnseq | paftools sam2paf - | pigz -p $cpu > ${pre}_reMEM.paf.gz
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> mem and paf: non-zero exit." >&2;exit $?;fi
	echo "generated MEM aln paf: ${pre}_reMEM.paf.gz" >&2
fi

echo "filtering raw unalign records out if coverage >= $tcvrg and identity >= $tpid ..." >&2
# check if 0 filter records
zcat ${pre}_reMEM.paf.gz |\
perl -lane '$F[0]=~s/\:\.//g;print $F[0] if $F[10]/$F[1] >= $ENV{tcvrg} && $F[9]/$F[10] >= $ENV{tpid}' |\
perl -lane 'print if ++$h{$_}==1' >${pre}_reMEM_filter.id  
num_reMEM_filter=$(cat ${pre}_reMEM_filter.id | wc -l)
echo "${pre}_reMEM_filter.id ---> $num_reMEM_filter " >&2
if [ $num_reMEM_filter -gt 0 ];then
	echo "extract filtered records to : ${pre}.unaligned.MEMfiltered.fa.gz" >&2
	seqkit replace -j $cpu -p ":. \d+" -r "" $unalnseq | seqkit grep -j $cpu -v -n -f ${pre}_reMEM_filter.id -o ${pre}.unaligned.MEMfiltered.fa.gz
else
	echo "No reMEM_filter records. Copy $unalnseq to ${pre}.unaligned.MEMfiltered.fa.gz" >&2
	seqkit replace -j $cpu -p ":. \d+" -r "" $unalnseq -o ${pre}.unaligned.MEMfiltered.fa.gz
fi


if [[ $? -ne 0 ]] ; then echo "[ERROR] --> filter reMEM: non-zero exit." >&2;exit $?;fi
echo "generated filtered unaligned seq: 
	${pre}.unaligned.MEMfiltered.fa.gz
	" >&2
# stats
seqkit stats -j $cpu $unalnseq ${pre}.unaligned.MEMfiltered.fa.gz >&2

