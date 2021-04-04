#!/usr/bin/env bash
usage="
------------------------------
cluster by bed and seq
------------------------------
USAGE:
	$0 <bed.gz> <unaln_seq_list> <step_pre> <distance> <coverage> <identity> <wordlength> <threads> <keep>
------------------------------
                   Songtao Gui
              songtaogui@sina.com
"
if [[ $# -ne 9 ]]; then 
	echo "$usage"
	exit 1
fi
export in_bed=$1
export unaln_seq_list=$2
export pre=$3
export pre_distance=$4
export coverage=$5
export identity=$6
export wordlength=$7
export threads=$8
export keep=$9

# --------------------
workdir=$PWD
srcdir=$(cd $(dirname $0); pwd)
cd $workdir;
export cdhit_script=$srcdir/pan11_run_cdhit_mtffa.sh
# ---------------------

#------------check and skip finished step-------
progress_check=$PWD/progress_check
#-----------------------------------------------
mkdir -p $progress_check
# 01
echo "[${pre}_STEP01]>>> merging records within $pre_distance bp ..." >&2
if [ ! -s $in_bed ];then
	echo "[ERROR] -->  no file $in_bed" >&2
	exit 1
fi
if [ -s "${progress_check}/${pre}_01.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_01.finished >&2
else
	sortBed -i $in_bed | mergeBed -i stdin -d $pre_distance -o count_distinct,distinct -c 4,4 | pigz -p $threads >${pre}_01_merged_${pre_distance}bp.bed.gz &&\
	ls ${pre}_01_merged_${pre_distance}bp.bed.gz >&2 
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 01 mergeBed: non-zero exit." >&2
		rm -f ${pre}_01_merged_${pre_distance}bp.bed.gz
		exit 1
	fi
	echo "Already finished ${pre}_01_merged_${pre_distance}bp.bed.gz" >${progress_check}/${pre}_01.finished
fi
# 02
echo "[${pre}_STEP02]>>> get un merged records ..." >&2 
if [ -s "${progress_check}/${pre}_02.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_02.finished >&2
else
	csvtk filter -j $threads -tTH -f "4<2" ${pre}_01_merged_${pre_distance}bp.bed.gz |\
	csvtk cut -tTH -f 5 -j $threads |\
	csvtk uniq -j $threads -o ${pre}_02_uniqID_uncls.id.gz &&\
	ls ${pre}_02_uniqID_uncls.id.gz >&2
	# BEref_SameOrd#5to3#L_R#CM009906.1_216895341_217264269_Mo17_scaftig0000608_264880-270505#++
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 02 Get uniqID: non-zero exit." >&2
		rm -f ${pre}_02_uniqID_uncls.id.gz
		exit 1
	fi
	echo "Already finished ${pre}_02_uniqID_uncls.id.gz" >${progress_check}/${pre}_02.finished
fi
# 03
# get clusterID of the rest cls
echo "[${pre}_STEP03]>>> Get clusterID to ${pre}_03_EACH_sample_CLS_tsv dir ..."
if [ -s "${progress_check}/${pre}_03.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_03.finished >&2
else
	csvtk filter -j $threads -tTH -f "4>1" ${pre}_01_merged_${pre_distance}bp.bed.gz | perl -lane '
	$,="\t";
	$cls=join("_","$ENV{pre}CLS",@F[0..3]);
	@id=split(/,/,$F[4]);
	foreach $id (@id){
		($sp,$tag,$seqid)=(split(/#/,$id))[0,1,4];
		if($sp && $tag && $seqid){
			print $cls,$sp,$id,$seqid
		}else{
			die "cannot parse features from $id";
			exit(1);
		}
	}' | csvtk split -j $threads -tTH -f 2 -o ${pre}_03_EACH_sample_CLS_tsv
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 03 get sample CLS : non-zero exit." >&2
		rm -f ${pre}_03_EACH_sample_CLS_tsv
		exit 1
	fi
	echo "Already finished ${pre}_03_EACH_sample_CLS_tsv" >${progress_check}/${pre}_03.finished
fi
# 04
echo "[${pre}_STEP04]>>> Get CLS sequences ... " >&2
if [ -s "${progress_check}/${pre}_04.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_04.finished >&2
else
	echo ">>>>>> get sequences by sample to ${pre}_04_EACH_sample_CLS_tfa dir ..." >&2
	mkdir -p ${pre}_04_EACH_sample_CLS_tfa 
	# get sequences by CLS
	######### run loop with multi threads #########
	# set threads
	loop_threads=$threads 
	# open fd 
	[ -e /tmp/fd1 ] || mkfifo /tmp/fd1  
	exec 7<>/tmp/fd1
	rm -rf /tmp/fd1
	for i in $(seq 0 $loop_threads);do
	    echo
	done >&7
	########### progress bar ###########
	####### scripts before loop ########
	export ttt=$(ls ${pre}_03_EACH_sample_CLS_tsv/*.tsv | wc -l)
	echo "Progress:( Total records: $ttt )
start: $(date)
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2
	# creat tmp files for counting (use tmp files instead of VARs for FIFO-loop)
	rm -f progress_bar_kkk.tmp progress_bar_ccc.tmp
	echo "start" > progress_bar_ccc.tmp
	touch progress_bar_kkk.tmp
	#####################################
	ls ${pre}_03_EACH_sample_CLS_tsv/*.tsv | while read cls_file
	do
	read -u7
	{	
		tpre=$(basename $cls_file)
		tpre=${tpre%%.tsv}
		tpre=${tpre##stdin-}
		cur_fa=$(grep -m 1 "/${tpre}.unaligned" $unaln_seq_list)
		if [ ! -s "$cur_fa" ];then
			echo "[ERROR] --> No file: $cur_fa" >&2
			exit 1
		fi
		csvtk join -tTH -f "4;1" -k $cls_file <(seqkit fx2tab $cur_fa) -o ${pre}_04_EACH_sample_CLS_tfa/${tpre}_cls_tfa.tsv.gz
		#### Progress bar within loop ####
		ccc=$(perl -e '$a=`cat progress_bar_ccc.tmp | wc -l`;print int($a*100/$ENV{ttt});')
		kkk=$(cat progress_bar_kkk.tmp | wc -l)
		echo "ccc" >> progress_bar_ccc.tmp
		if [[ $ccc -eq $kkk ]];then
			echo -e "kkk1\nkkk2" >> progress_bar_kkk.tmp
			echo -n "*" >&2;
		fi
		###################################
		echo >&7
	} &
	done | cat
	wait
	exec 7>&- # close fd
	exec 7<&- # close fd
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 04 seq tfa: non-zero exit."
		exit 1
	fi
	echo -e "\nDone." >&2
	echo ">>>>>> Get CLS sequences by Cluster to ${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz ..." >&2

	cat  ${pre}_04_EACH_sample_CLS_tfa/*_cls_tfa.tsv.gz |pigz -d -p $threads |\
	perl -lane '$a=join("__TAB__",@F[1..$#F]);print "$F[0]\t$a";'|\
	csvtk sort -tTH -k 1 -j $threads | csvtk collapse -tTH -j $threads -f 1 -v 2 -s "," -o ${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 04 cls tfa: non-zero exit." >&2
		exit 1
	fi	
	echo "Done." >&2
	ls -lh ${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz >&2

	# cat ${pre}_04_EACH_sample_CLS_tfa/*_cls_tfa.tsv.gz | csvtk split -tTH -j $threads -f 1 -G -o ${pre}_04_EACH_cluster_CLS_tfa
	# if [[ $? -ne 0 ]] ; then 
	# 	echo "[ERROR] --> 04 cls tfa: non-zero exit." >&2
	# 	exit 1
	# fi

	# remove intermediate files if $keep is not "T"
	if [ "$keep" != "T" ];then
		echo "cleaning ${pre}_04_EACH_sample_CLS_tfa dir ... " >&2
		rm -rf ${pre}_04_EACH_sample_CLS_tfa
		echo "Done." >&2
	fi
	echo "Already finished ${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz" >${progress_check}/${pre}_04.finished
fi

# 05
# run cdhit for each CLS 
workdir=$PWD
if [ -s "${workdir}/${pre}_05_All_cdhit_sequences.fa.gz" -a -s "${workdir}/${pre}_05_All_cdhit_clstr.tsv.gz" ];then
	echo "Already finished ${pre}_05_All_cdhit_clstr.tsv.gz and ${pre}_05_All_cdhit_sequences.fa.gz" >${progress_check}/${pre}_05.finished
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_05.finished >&2
else
	############# USE LSF or PBS ############
	echo "[${pre}_STEP05]>>> Split${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz into $threads parts and run each using LSF or PBS ... " >&2
	echo "Output part files to ${pre}_05_split_run ..." >&2
	mkdir -p ${pre}_05_split_run
	cd ${pre}_05_split_run
	# get total number
	if [ -s "ttt.num" ];then
		export ttt=$(cat ttt.num)
	else
		pigz -p $threads -dc $workdir/${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz | wc -l >ttt.num
		if [[ $? -ne 0 ]] ; then 
			echo "[ERROR] --> ttt.num: non-zero exit." >&2
			rm -f ttt.num
			exit 1
		fi
		export ttt=$(cat ttt.num)
	fi
	split_num=$((ttt/threads+1))
	# get split files
	echo "split with -N $split_num ... " >&2
	pigz -dc -p $threads $workdir/${pre}_04_EACH_cluster_CLS_mttfa.tsv.gz | parallel -j $threads --pipe -N $split_num 'cat > split_mttfa_{#}.tsv' 
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> split-mttfa: non-zero exit." >&2
		exit 1
	fi
	echo "Done!" >&2
	# echo "Please run the script below for each splitted mttfa.tsv (using LSF or PBS if you like):
	# >>> run_cdhit_mtffa.sh <mttfa.tsv> ${pre} $coverage $identity $wordlength $keep
	# 	and merge the outputed tsv and fa to:
	# >>> ${workdir}/${pre}_05_All_cdhit_clstr.tsv.gz
	# >>> ${workdir}/${pre}_05_All_cdhit_sequences.fa.gz
	# 	then rerun this script to continue.
	# " >&2
	# exit;
	# USAGE:0 input1 input2
	cdhit_mtffa () {
		# NOTE: run cdhit for each split file with parallel
		local mttfa=$1
		bash $cdhit_script ${mttfa} ${pre} ${coverage} ${identity} ${wordlength} ${keep} 1> ${mttfa}.log 2>&1
	}
	export -f cdhit_mtffa
	# DEBUG: parallel run cdhit_mtffa -> test if work
	parallel -j $threads --bar cdhit_mtffa :::: <(ls split_mttfa_*.tsv)
	if [ $? -ne 0 ];then gst_err "parallel run cdhit mtffa failed: Non-zero exit"; exit 1;fi
	# NOTE: merge split result to concated output and clean intimate files
	cat *.clstr | pigz -p $threads > ${workdir}/${pre}_05_All_cdhit_clstr.tsv.gz &&\
	cat *.fa | pigz -p $threads > ${workdir}/${pre}_05_All_cdhit_sequences.fa.gz
	if [ $? -ne 0 ];then gst_err "get merged results failed: Non-zero exit"; rm -f ${workdir}/${pre}_05_All_cdhit_clstr.tsv.gz ${workdir}/${pre}_05_All_cdhit_sequences.fa.gz; exit 1;fi
	# NOTE: clean split dir when every thing finished
	cd ../ && rm -rf ${pre}_05_split_run
fi

# 06
echo "[${pre}_STEP06]>>> collapse cdhit clstr records ..." >&2
cd ${workdir}
if [ -s "${progress_check}/${pre}_06.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_06.finished >&2
else
	#csvtk collapse -tTH -j $threads -f 1 -s "," -v 2 ${pre}_05_All_cdhit_clstr.tsv.gz -o ${pre}_06_collapsed_cdhit_clstr.tsv.gz
	zcat ${pre}_05_All_cdhit_clstr.tsv.gz | bedtools groupby -g 1 -c 2,2,2 -o first,count_distinct,distinct | pigz -p $threads > ${pre}_06_collapsed_cdhit_clstr.tsv.gz 
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 06 collapse clstr: non-zero exit." >&2
		rm -f ${pre}_06_collapsed_cdhit_clstr.tsv.gz
		exit 1
	fi
	echo "Done. File in ${pre}_06_collapsed_cdhit_clstr.tsv.gz
	Format: CLS_ID | Ref_seq_ID | Number_Seqs | All_Seqs(comma)" >&2
	echo "Already finished ${pre}_06_collapsed_cdhit_clstr.tsv.gz" >${progress_check}/${pre}_06.finished
fi

# 07
# get ref 
echo "[${pre}_STEP07]>>> Get ref record ID of each clstr ..." >&2
if [ -s "${progress_check}/${pre}_07.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/${pre}_07.finished >&2
else
	csvtk grep -tTH -j $threads -f 3 -p "00ref00" ${pre}_05_All_cdhit_clstr.tsv.gz |\
	csvtk cut -tTH -f 2 -j $threads -o ${pre}_07_clstr_refID.id.gz 
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> 07 get ref id: non-zero exit." >&2
		rm -f ${pre}_07_clstr_refID.id.gz 
		exit 1
	fi
	echo "Done. File in ${pre}_07_clstr_refID.id.gz " >&2
	echo "Already finished ${pre}_07_clstr_refID.id.gz" >${progress_check}/${pre}_07.finished
fi

rm -f progress_bar_ccc.tmp progress_bar_kkk.tmp