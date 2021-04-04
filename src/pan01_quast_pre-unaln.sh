#!/usr/bin/env bash
if [ $# -ne 6 ]; then echo "usage $0 <pre> <scaftig.fa> <ref.fa> <minlen> <scftig_len> <cpu>" >&2; exit 1;fi

pre=$1
scft_fa=$2
ref=$3
export minlen=$4
scftig_len=$5
cpu=$6

echo "
>>>> run QUAST5 
" >&2
if [ ! -s $scft_fa ];then echo "[ERROR]---> No scaftig file: $scft_fa" >&2;exit 1;fi

unalninfo=${pre}_quast/contigs_reports/contigs_report_${pre}_scaftig${scftig_len}.unaligned.info
fcoords=${pre}_quast/contigs_reports/minimap_output/${pre}_scaftig${scftig_len}.coords.filtered

if [ -s "$unalninfo" -a -s "$fcoords" ];then 
	echo "Skip running QUAST5. Output files exist: ${pre}_quast " >&2
else
	echo "Running Quast 5 ..." >&2
	quast.py $scft_fa --debug --eukaryote --fragmented --no-snps --no-plots --no-icarus --no-sv -R $ref --threads $cpu -o ${pre}_quast 
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> quast: non-zero exit." >&2;exit $?;fi
    echo "Quast5 finished nomarlly." >&2
fi

echo "
>>>> get ump seq from quast unalign info " >&2
if [ -s "${pre}.unaligned.${minlen}bp.fa.gz" -a -s "${pre}.unaligned.${minlen}bp.bed" ];then
	echo "Skip running. Output files exist: 
	${pre}.unaligned.${minlen}bp.bed
	${pre}.unaligned.${minlen}bp.fa.gz
	" >&2;	
else
	# get_unaligned_fasta.sh $scft_fa $unalninfo $minlen ${pre}
	# get > ${minlen} bp bed (0 based):
	if [ -s $unalninfo ];then
		## get unalign percentage and unalign fragment length
	 	echo "get unaligned bed > ${minlen} bp ..." >&2
		sed '1d' $unalninfo| perl -lane '@un=split(/,/,$F[-1]);$,="\t";foreach $uu (@un){ ($s,$e)=split(/-/,$uu);$len=$e-$s+1; print $F[0],$s-1,$e,"unaln",$len,"+" if $len >= $ENV{minlen};}' >${pre}.unaligned.${minlen}bp.bed
	    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get unaln bed: non-zero exit." >&2;exit $?;fi
	else
		echo "[ERROR]---> No $unalninfo !" >&2;exit;
	fi
	# get seq by bed:
	if [ -s $scft_fa ];then
		echo "get unaligned fasta > ${minlen} bp ..." >&2
		seqkit subseq -j $cpu --bed ${pre}.unaligned.${minlen}bp.bed $scft_fa -o - | seqkit replace -j $cpu -p ":.*" -r "" -o ${pre}.unaligned.${minlen}bp.fa.gz
		# seqkit subseq -j $cpu --bed ${pre}.unaligned.${minlen}bp.bed $scft_fa -o - | seqkit fx2tab -j $cpu | perl -F"\t" -lane '$F[0]=~s/:.*//;$,="\t";print @F;' | seqkit tab2fx -j $cpu -o ${pre}.unaligned.${minlen}bp.fa.gz
	    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get unaln seq: non-zero exit." >&2;exit $?;fi
	    echo "generated unaligned fasta: ${pre}.unaligned.${minlen}bp.fa.gz " >&2
	else
		echo "[ERROR]---> No $scft_fa  !" >&2;exit;
	fi
fi
# stats
wc -l ${pre}.unaligned.${minlen}bp.bed >&2
seqkit stats -j $cpu ${pre}.unaligned.${minlen}bp.fa.gz >&2