#!/usr/bin/env bash
usage="
------------------------------
RUN cdhit in one threads for each split mttfa.
Use this script with PBS or LSF system to avoid out of memory.
Used for PANZ_cluster only.
------------------------------
USAGE:
	$0 <mttfa> <pre:S1_S2_S3> <coverage> <identity> <wordlength> <Keep_tmp: T or F>
------------------------------
                   Songtao Gui
              songtaogui@sina.com
"
if [[ $# -ne 6 ]]; then 
	echo "$usage"
	exit 1
fi

mttfa_file=$1
pre=$2
threads=1
coverage=$3
identity=$4
wordlength=$5
keep=$6
workdir=$PWD
rm -rf ${mttfa_file%%.tsv}_cdhit_out
mkdir ${mttfa_file%%.tsv}_cdhit_out
cd ${mttfa_file%%.tsv}_cdhit_out
if [ ! -s "../$mttfa_file" ];then
	echo "[ERROR] --> No file: $PWD/../$mttfa_file" >&2
	exit 1
fi
ttt=$(cat $workdir/$mttfa_file | wc -l)
echo "
CDHIT CMD: cd-hit-est -d 0 -aS $coverage -c $identity -n $wordlength -G 0
Progress:( Total records: $ttt )
start: $(date)
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2

# creat tmp files for counting (use tmp files instead of VARs for FIFO-loop)
rm -f progress_bar_kkk.tmp progress_bar_ccc.tmp
echo "start" > progress_bar_ccc.tmp
touch progress_bar_kkk.tmp
#####################################
# ######### run loop with multi threads #########
# # set threads
# loop_threads=$threads
# # open fd 
# [ -e /tmp/fd1 ] || mkfifo /tmp/fd1  
# exec 7<>/tmp/fd1
# rm -rf /tmp/fd1
# for i in $(seq 0 $loop_threads);do
#     echo
# done >&7
# ##################################
####### tmp file name ##########
tmp_file=Combine_out_$(basename ${mttfa_file%%.tsv})
#echo $tmp_file
################################
rm -f TEMP_out_*
cat $workdir/$mttfa_file | while read tpre mttfa
do
# read -u7
# {	
	#### Progress bar within loop ####
	line_num=$(cat progress_bar_ccc.tmp | wc -l)
	ccc=$((line_num*100/ttt))
	kkk=$(cat progress_bar_kkk.tmp | wc -l)
	echo "ccc" >> progress_bar_ccc.tmp
	if [[ $ccc -eq $kkk ]];then
		echo -e "kkk1\nkkk2" >> progress_bar_kkk.tmp
		echo -n "*" >&2;
	fi
	echo -e "$tpre\t$mttfa" | perl -lane 's/,/\n$F[0]\t/g;s/__TAB__/\t/g;print ;'|\
	perl -lane 'print "$F[0]\@$F[2]\t$F[4]" if $#F==4;
		#die "Wrong Number of column for $F[2]" if $#F!=4;
	' | seqkit tab2fx >${tpre}.tmp.fa
	if [ ! -s "${tpre}.tmp.fa" ];then
		echo "[ERROR] --> No file: ${tpre}.tmp.fa" >&2
		#exit 1
	fi
	echo "cd-hit-est -aS $coverage -c $identity -n $wordlength -d 0 -G 0 -i ${tpre}.tmp.fa -o ${tpre}_cdhit.fa" >>${tmp_file}.cdhit_logs
	cd-hit-est -aS $coverage -c $identity -n $wordlength -d 0 -G 0 -i ${tpre}.tmp.fa -o ${tpre}_cdhit.fa 1>/dev/null 2>&1 &&\
	cur_cls_tsv=$( cat ${tpre}_cdhit.fa.clstr | perl -F">" -lane '
		if(/^>Cluster (\d+)/){
			$cls="cdhit_$1";
		}else{
			($lc,$id)=split(/@/,$F[-1]);
			$id=~s/\.\.\.\sat\s/\t/g;
			$id=~s/\.\.\.\s\*/\t00ref00/g;
			if($id){
				printf "%s_%s\t%s\n",$lc,$cls,$id;
			}else{
				die "[ERROR] --> No id for $_";
				exit(1);
			}
		} ' | sort -k1,1 -k3,3 ) &&\
	cur_cls_fa=$( cat ${tpre}_cdhit.fa ) &&\
	echo "$cur_cls_tsv" >> ${tmp_file}.tsv &&\
	cat ${tpre}_cdhit.fa.clstr >> ${tmp_file}.clstr &&\
	cat ${tpre}_cdhit.fa >>  ${tmp_file}.fa 
		# format: | CLS_ID | seqID | cls_aln_info(ref or identity) |
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> cdhit: non-zero exit: $tpre" >&2
		exit 1
	fi
	# clean each cdhit file unless $keep == T
	if [ "$keep" != "T" ];then
		rm -f ${tpre}.tmp.fa ${tpre}_cdhit.fa ${tpre}_cdhit.fa.clstr
	fi		
	# echo  >&7
# } &
done | cat
wait 
# exec 7>&- # close fd
# exec 7<&- # close fd
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> cdhit: non-zero exit." >&2
	exit 1
fi
echo "done! Out put to $(ls ${tmp_file}*)"
