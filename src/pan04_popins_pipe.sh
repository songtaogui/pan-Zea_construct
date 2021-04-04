#!/usr/bin/env bash
# run popins for each sample

if [ $# -ne 7 ]; then echo "Usage: $0 <prefix> <raw_bam> <reMEM_filter_scftg> <reference> <read_length> <cpu> <minlen> " >&2; exit 1;fi

pre=$1
bam=$2
ctggz=$3
ref=$4
rlen=$5
cpu=$6
minlen=$7
mem=2G
ctg_pre=clean_nonrefseq_${minlen}_$pre
ctg=${ctg_pre}.fa
dir=.

if [ "${ctggz:0:1}" != "/" ];then 
	echo "[ERROR] --> [-1]: File should use ABSOLUTE PATH: $ctggz" >&2;
	exit 1;
fi


if [ ! -s $ctggz ];then echo "[ERROR] --> No $ctggz found." >&2;exit 1;fi
if [ ! -s $ref ];then echo "[ERROR] --> No $ref found." >&2;exit 1;fi

mkdir -p pp_${pre}/${pre};
cd pp_${pre};

echo "# ---> Working dir: $PWD" >&2

if [ -s $bam ];then echo "Loaded $bam ." >&2;else echo "[ERROR] --> No $bam found." >&2;exit 1;fi
if [ ! -s ${bam}.bai ]; then 
	echo "Index $bam ." >&2
	samtools index $bam
fi


echo "[`date +"%m-%d %H:%M"`] ----------> popins get unmapped reads" >&2
if [[ -s "${pre}/non_ref.bam" ]];then
	echo "skip get unmapped reads. File exists." >&2
else
	if [ -s $bam ];then echo "Loaded $bam ." >&2;else echo "[ERROR] --> No $bam found." >&2;exit 1;fi
	popins assemble -p $dir -s $pre -t $cpu -m $mem $bam 
	# no velvet, would pop an ERROR report. ignore it.
	#if [[ $? -ne 0 ]] ; then echo "[ERROR] --> popins assemble : non-zero exit." >&2;exit $?;fi
fi

if [ ! -s $ctg ];then
	echo "Get ctgfile: $ctg from $ctggz" >&2
	pigz -p $cpu -d -c $ctggz > $ctg
	if [[ $? -ne 0 ]] ; then echo "## get ctg --> returned non-zero exit status - aborting" >&2;exit $?;fi
fi

if [ -d "$pre" ];then
	echo "Loaded $pre dir." >&2
else
	echo "[ERROR] --> No $pre dir found." >&2
	exit 1;
fi

if [ -s "$ref" -a -s "$ctg" ];then
	echo "Loaded $ref and $ctg ." >&2
else
	echo "[ERROR] --> No $ref or $ctg !" >&2
	exit 1:
fi

echo "[`date +"%m-%d %H:%M"`] ----------> popins contigmap" >&2
if [ -s "${pre}/non_ref_new.bam" -a -s "${pre}/non_ref_new.bam.bai" ];then
	echo "SKIP. File exists: $PWD/${pre}/non_ref_new.bam" >&2
else
	popins contigmap -p $dir -c $ctg -r $ref -t $cpu -m $mem --best --maxInsertSize 800 $pre >&2
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> popins contigmap returned non-zero exit status - aborting" >&2;exit 1;fi
fi

if [ ! -s "${ctg_pre}_genotype.FinishedLog" ];then
	echo "[`date +"%m-%d %H:%M"`] ----------> popins place-refalign"	 >&2
	popins place-refalign -p $dir -c $ctg -r $ref -l ${ctg_pre}_locations.txt -i ${ctg_pre}_insertions.vcf -g ${ctg_pre}_groups.txt --minScore 0.3 --minReads 2 --maxInsertSize 800 --readLength $rlen --groupDist 100
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> popins place-refalign returned non-zero exit status - aborting" >&2;exit $?;fi

	echo "[`date +"%m-%d %H:%M"`] ----------> popins place-splitalign" >&2		
	popins place-splitalign -p $dir -c $ctg -r $ref --maxInsertSize 800 --readLength $rlen $pre
	if [[ $? -ne 0 ]] ; then echo "[ERROR] -->  popins place-splitalign returned non-zero exit status - aborting" >&2;exit $?;fi

	echo "[`date +"%m-%d %H:%M"`] ----------> popins place-finish" >&2		
	popins place-finish -p $dir -r $ref -i ${ctg_pre}_insertions.vcf
	if [[ $? -ne 0 ]] ; then echo "[ERROR] -->  popins place-finish returned non-zero exit status - aborting" >&2;exit $?;fi

	echo "[`date +"%m-%d %H:%M"`] ----------> popins genotype" >&2
	#nf=`grep -m1 "#CHROM" ${ctg_pre}_insertions.vcf | awk '{print NF}'`

	popins genotype -i ${ctg_pre}_insertions.vcf -p $dir -c $ctg -r $ref -m RANDOM -w 50 $pre 
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> popins genotype returned non-zero exit status - aborting" >&2;exit $?;fi
	
	if [ -s ${pre}/insertions.vcf ];then
		perl -lane '/^#/ && print && next;$,="\t";print @F[0..8,10];' ${pre}/insertions.vcf > ${pre}/${ctg_pre}_genotyped_insertions.vcf
	else
		exit 1;
	fi
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> fmt genovcf returned non-zero exit status - aborting" >&2;exit $?;fi

	echo "Finished genotype for ${ctg_pre}_insertions_edited.vcf" >${ctg_pre}_genotype.FinishedLog
else
	echo "Skip genotype for ${ctg_pre}_insertions_edited.vcf. File exists." >&2
fi

echo "[`date +"%m-%d %H:%M"`] ----------> Done !" >&2



