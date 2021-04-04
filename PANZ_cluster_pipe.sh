#!/usr/bin/env bash 

usage="
--------------------------------------------------
 Clustering PANZ_individual_pipe.sh outputs
--------------------------------------------------
USAGE: 
    bash $0 <OPTIONS>
OPTIONS: (R:required  O:optional)
  -h                  show help and exit.
  -t <num>     O      set threads (default 2)
  -V <list>    R      VCF file list
  -U <list>    R      UNANCHORED list
  -S <list>    R      unaln sequences
  -d <num>     O      pre distance (default 200 bp)
  -D <num>     O      final distance (default 2000 bp)
  -c <0-1>     O      coverage cutoff (default 0.8)
  -i <0-1>     O      identity cutoff (default 0.9)
  -w <num>     O      cdhit word length (default 10)
  -n <num>     O      Min support number (default 3)
  -r <num>     O      Min REF support number (default 1)
  -k <T/F>     O      Keep all intermediate files if 'T'(default F)
--------------------------------------------------
                                    Songtao Gui
                                songtaogui@sina.com
"

############### Set Default OPT #############
threads=2
export VCF_list=
export UNANCHORED_list=
export unaln_seq_list=
pre_distance=200
final_distance=2000
export coverage=0.8
export identity=0.9
wordlength=10
export supp_num=3
export ref_num=1
keep=F
# --------------------
workdir=$PWD
srcdir=$(cd $(dirname $0); pwd)
cd $workdir;
# ---------------------
################# GETOPT ####################
if [ $# -eq 0 ]; then echo "$usage" >&2; exit 1;fi
while getopts t:V:U:S:d:D:c:i:w:n:r:k:h opt; 
do 
case $opt in
	h)  echo "$usage" >&2
		exit
		;;
	t)  threads=$OPTARG
		;;
	V)  VCF_list=$OPTARG
		;;
	U)  UNANCHORED_list=$OPTARG
		;;
	S)  unaln_seq_list=$OPTARG
		;;
	d)  pre_distance=$OPTARG
		;;
	D)  final_distance=$OPTARG
		;;
	c)  coverage=$OPTARG
		;;
	i)  identity=$OPTARG
		;;
	w)  wordlength=$OPTARG
		;;
	n)  supp_num=$OPTARG
		;;
	r)  ref_num=$OPTARG
		;;				
	k)  keep=$OPTARG
		;;
    :)  printf "[ERROR] --> missing argument for -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
    \?) printf "[ERROR] --> illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
esac
done
shift $((OPTIND - 1))
################# GETOPT END ################

echo "
#################  Pre-step  ######################
                generate bed files 
###################################################
" >&2
if [[ -s "00_combine-left.bed.gz" && -s "00_combine-right.bed.gz" ]];then
	echo -e "Skip. File exists:" >&2
else
	echo ">>> get left and right bed files from VCF file to 00_bed_files ..." >&2
	mkdir -p 00_bed_files
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

	cat $VCF_list | while read vcf_file
	do
	read -u7
	{
		if [ ! -s $vcf_file ];then
			echo "No file: $vcf_file" >&2
			exit 1
		fi
		# get basename
		export vcf_pre=$(basename $vcf_file)
		vcf_pre=${vcf_pre%%_00_FINAL_ALL.vcf}
		cat $vcf_file | perl -F"\t" -lane '
			$,="\t";
			/^#/ && next;
			($chr1,$end1,$id)=@F[0,1,2];
			$start1=$end1-1;
			$start1=0 if $start1 < 0;
			if(/CHR2=(.*);END=(\d+);/){
				$chr2=$1;
				$end2=$2;
				$start2=$end2-1;
				$start2=0 if $start2 < 0;
				print $chr1,$start1,$end1,"$ENV{vcf_pre}#$id","left";
				print $chr2,$start2,$end2,"$ENV{vcf_pre}#$id","right";
			}else{
				die "No CHR2 or END for $id";
			}
		' >00_bed_files/${vcf_pre}.bed
		if [[ $? -ne 0 ]] ; then 
			echo "[ERROR] --> get-bed: non-zero exit." >&2
			rm -f 00_bed_files/${vcf_pre}.bed
			exit 1
		fi
		echo >&7
	} &
	done | cat 
	wait 
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> get-bed: non-zero exit." >&2
		exit 1
	fi
	exec 7>&- # close fd
	exec 7<&- # close fd

	echo ">>> generate left and right bed: 00_combine-left.bed.gz 00_combine-right.bed.gz ..." >&2
	pwd

	csvtk concat -j $threads -tTH 00_bed_files/*.bed | csvtk split -G -j $threads -tTH -f 5 -o $PWD
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> csvtk split: non-zero exit." >&2
		# rm -f xxx.outputs
		exit 1
	fi
	echo "Done csvtk split" >&2
	mv stdin-left.tsv.gz 00_combine-left.bed.gz &&\
	mv stdin-right.tsv.gz 00_combine-right.bed.gz &&\
	echo "Done. Generated bed files:" >&2
fi
ls -lh 00_combine-left.bed.gz 00_combine-right.bed.gz >&2

echo "
###################  STEP1  #######################
 Clustering by left-posi within $pre_distance bp 
###################################################
" >&2
echo "bash $srcdir/src/pan10_bed_cls.sh 00_combine-left.bed.gz $unaln_seq_list S1 $pre_distance $coverage $identity $wordlength $threads $keep" >&2
bash $srcdir/src/pan10_bed_cls.sh 00_combine-left.bed.gz $unaln_seq_list S1 $pre_distance $coverage $identity $wordlength $threads $keep

if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> step1: non-zero exit."
	exit 1
fi
echo "Done. Output files:" >&2
ls | grep "^S1" >&2
# bash panz_bed_cls.sh <bed> <unaln_seq_list> <step_pre> <distance> <coverage> <identity> <wordlength> <threads> <keep>

echo "
###################  step2   ######################
 Clustering by right-posi within $pre_distance bp 
###################################################
"
# get right cls beds
if [[ ! -s "S1_02_uniqID_uncls.id.gz" || ! -s "S1_07_clstr_refID.id.gz" ]];then
	echo "No file: S1_02_uniqID_uncls.id.gz OR S1_07_clstr_refID.id.gz" >&2
	exit 1
fi
echo ">>> generate bed file for step2: S2_00_combine-right.bed.gz ..."
if [ -s "S2_00_combine-right.bed.gz" ];then
	echo "Skip. File exists." >&2
else
	cat S1_02_uniqID_uncls.id.gz S1_07_clstr_refID.id.gz | csvtk uniq -tTH -j $threads -o S1_10_all_ID_for_S2.id.gz &&\
	csvtk join -tTH -j $threads -f '4;1' 00_combine-right.bed.gz S1_10_all_ID_for_S2.id.gz -o S2_00_combine-right.bed.gz 
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> generate bed: non-zero exit." >&2
		rm -f S2_00_combine-right.bed.gz
		exit 1
	fi
	echo "Done. Generated file: S2_00_combine-right.bed.gz" >&2
fi

bash $srcdir/src/pan10_bed_cls.sh S2_00_combine-right.bed.gz $unaln_seq_list S2 $pre_distance $coverage $identity $wordlength $threads $keep
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> step2: non-zero exit."
	exit 1
fi
echo "Done. Output files:" >&2
ls | grep "^S2" >&2

 # final cluster
echo "
###################  step3   ######################
   Final clustering within $final_distance 
###################################################
"
# get final cls beds
if [[ ! -s "S2_02_uniqID_uncls.id.gz" || ! -s "S2_07_clstr_refID.id.gz" ]];then
	echo "No file: S2_02_uniqID_uncls.id.gz OR S2_07_clstr_refID.id.gz" >&2
	exit 1
fi
echo ">>> generate bed files for step3: " >&2
if [ -s "S3_00_all_left_right.bed.gz" ];then
	echo "Skip. File exists." >&2
else
	cat S2_02_uniqID_uncls.id.gz S2_07_clstr_refID.id.gz | csvtk uniq -tTH -j $threads -o S2_10_all_ID_for_S3.id.gz &&\
	cat 00_combine-left.bed.gz 00_combine-right.bed.gz > 00_all_left_right.bed.gz &&\
	csvtk join -tTH -j $threads -f '4;1' 00_all_left_right.bed.gz S2_10_all_ID_for_S3.id.gz -o S3_00_all_left_right.bed.gz
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> generate bed: non-zero exit." >&2
		rm -f S3_00_all_left_right.bed.gz
		exit 1
	fi
	echo "Done. Generated file: S3_00_all_left_right.bed.gz" >&2
fi

bash $srcdir/src/pan10_bed_cls.sh S3_00_all_left_right.bed.gz $unaln_seq_list S3 $final_distance $coverage $identity $wordlength $threads $keep
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> step3: non-zero exit."
	exit 1
fi
echo "Done. Output files:" >&2
ls | grep "^S3" >&2

echo "
################ STEP4 ############################
 Merge all clusters and filter
###################################################
"
# merge S*_06_collapsed_cdhit_clstr.tsv.gz
echo "[S4_STEP01]>>> merge 3 cluster results and get final cluster file ..." >&2
if [ -s S4_00_FINAL_Cluster.list.gz ];then 
	echo "Skip. File exists" >&2
	else
	############ Check if all related files exists ###########
	num_related_file=1;
	for related_file in S*_06_collapsed_cdhit_clstr.tsv.gz
	do
		if [[ ! -s "$related_file" ]]; then
			printf "[ERROR] --> No file: %s \n" $related_file >&2
			let num_related_file++
		fi
	done
	if [ $num_related_file -eq 1 ];then
		echo "All related files were found. Proceeding ..." >&2
	else
		echo "Check if you miss somthing." >&2
		exit 1;
	fi
	############### merging cluster files #################
	zcat S3_06_collapsed_cdhit_clstr.tsv.gz | cut -f 4 | perl -lane '
	BEGIN{
		$fs="S2_06_collapsed_cdhit_clstr.tsv.gz";
		open(IN,"gzip -dc $fs |") or die("cannot open $fs:$!");
		while(<IN>){
			chomp;
			(undef,$ref,undef,$all)=split(/\t/,$_);
			die("No record for $_") unless $ref;
			$hs{$ref}=$all;
		}
		close IN;
	}
	@seqid=split(/,/,$_);
	for ($i=0;$i<=$#seqid;$i++) {
		$seqid[$i]=$hs{$seqid[$i]} if $hs{$seqid[$i]};
	}
	$,=",";
	print @seqid;
	' | perl -lane '
	BEGIN{
		$fs="S1_06_collapsed_cdhit_clstr.tsv.gz";
		open(IN,"gzip -dc $fs |") or die("cannot open $fs:$!");
		while(<IN>){
			chomp;
			(undef,$ref,undef,$all)=split(/\t/,$_);
			die("No record for $_") unless $ref;
			$hs{$ref}=$all;
		}
		close IN;
	}
	@seqid=split(/,/,$_);
	for ($i=0;$i<=$#seqid;$i++) {
		$seqid[$i]=$hs{$seqid[$i]} if $hs{$seqid[$i]};
	}
	$,=",";
	print @seqid;
	' | perl -F"," -lane '
	%hash={};
	@uniq=grep {++$hash{$_}==1} @F;
	$cls=join(",",@uniq);
	$,="\t";
	print $#uniq+1,$cls;
	' | pigz -p $threads > S4_00_FINAL_Cluster.list.gz
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> merging cls: non-zero exit." >&2
		rm -f S4_00_FINAL_Cluster.list.gz
		exit 1
	fi
fi
echo "Done. Generated merged cluster file: S4_00_FINAL_Cluster.list.gz"
############## Filtering merged clusters ##########
# r>0  
#    t > r -> PASS_WithRef
#    t = r -> PASS_OnlyRef
#    t < r -> error
# r=0 
#    t = nlowqual -> LOWQUAL
#    t > nlowqual -> PASS_NoRef
#    t < nlowqual -> error
# r<0  error
Filter_CLS_fun(){	
	cur_tsupp=$(echo "$1" | cut -f 1 )
	cur_rcd=$(echo "$1" | cut -f 2 )
	cur_rsupp=$(echo "$cur_rcd" | grep -o "ref" | grep -c "ref" )
	nlowqual=$(echo "$cur_rcd" | grep -o "LOWQUAL" | grep -c "LOWQUAL" )
	if [[ "$cur_tsupp" -ge "$supp_num" || "$cur_rsupp" -ge "$ref_num" ]];then
		if [[ "$cur_rsupp" -gt 0 ]];then
			# r > 0	
			if [[ "$cur_tsupp" -gt "$cur_rsupp" ]];then
				#    t > r -> PASS_WithRef
				cls_type="PASS_WithRef"
			elif [[ "$cur_tsupp" -eq "$cur_rsupp" ]]; then
				#    t = r -> PASS_OnlyRef
				cls_type="PASS_OnlyRef"
			else
				#    t < r -> error
				echo -e "[ERROR] --> ref_num gt total num:\n $cur_tsupp\t$cur_rcd" >&2
				exit 1
			fi
		elif [[ "$cur_rsupp" -eq 0 ]];then
			# r = 0
			if [[ "$cur_tsupp" -gt "$nlowqual" ]];then
				#    t > nlowqual -> PASS_NoRef
				cls_type="PASS_NoRef"
			elif [[ "$cur_tsupp" -eq "$nlowqual" ]]; then
				#    t = nlowqual -> LOWQUAL
				cls_type="LOWQUAL"
			else
				#    t < nlowqual -> error
				echo -e "[ERROR] --> nlowqual gt total num:\n $cur_tsupp\t$cur_rcd" >&2
				exit 1
			fi
		elif [[ "$cur_rsupp" -lt 0 ]]; then
			# r < 0
			echo -e "[ERROR] --> ref num less than zero:\n $cur_tsupp\t$cur_rcd" >&2
			exit 1
		fi
		cur_out_line=$(echo -e "$cls_type\t$cur_tsupp\t$cur_rsupp\t$nlowqual\t$cur_rcd") &&\
		echo "$cur_out_line"
	fi
}
export -f Filter_CLS_fun

echo "[S4_STEP02]>>> Filtering out clusters with less than $supp_num total_supports or less than $ref_num Ref_supports ..." >&2
if [ -s S4_01_Filter_Cluster_fmt.tsv.gz ];then
	echo "SKip. File exists." >&2
else
	rm -f S4_01_Filter_Cluster_fmt.tsv
	pigz -dc -p $threads S4_00_FINAL_Cluster.list.gz | parallel -j $threads -k --lb --bar Filter_CLS_fun | pigz -p $threads >S4_01_Filter_Cluster_fmt.tsv.gz
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> filter by parallel: non-zero exit." >&2
		rm -f S4_01_Filter_Cluster_fmt.tsv.gz
		exit 1
	fi
fi
echo -e "Done. Filtered file:\n$(ls -lh S4_01_Filter_Cluster_fmt.tsv.gz)" >&2

echo "
[S4_STEP02]>>> Get represent sequences for each cluster ...
(Criteria: Length > BE|SE > REF|NonREF > PASS|LOWQUAL )
" >&2

# function
# input is S4_01_Filter_Cluster_fmt.tsv.gz per line
get_rep_cls_fun(){
	line=$1
	cls_type=$(echo "$line" | cut -f 1 )
	fmt_rcd=$(
		echo "$line" | cut -f 5 | sed 's/,/\n/g' | perl -lane '
			$line=$_;
			($sample,$type,$orient,$l_r,$id,$strand)=split(/#/,$line);
			$,="\t";
			if($id=~/_(\d+)-(\d+)$/){
				$length=$2-$1+1;
				$bs="SE";
				$bs="BE" if $type=~/^BE/;
				$ref=2;
				$ref=1 if $type=~/ref/;
				$pass=1;
				$pass=2 if $type=~/LOWQUAL/;
				print $length,$bs,$ref,$pass,$line;
			}else{
				die("Wrong scaffold id format: $id \n");
			}
		' )

		if [[ -z "$fmt_rcd" ]];then
			echo "[ERROR]-->empty fmt_rcd." >&2
			#exit 1
		fi

		rep_rcd=$(echo "$fmt_rcd" | sort -k1,1nr -k2,2 -k3,3n -k4,4n | cut -f 5 | head -1)
		
		if [[ -z "$rep_rcd" ]];then
			echo "[ERROR]-->empty rep_rcd." >&2
			#exit 1
		fi

		echo -e "$rep_rcd\t$line"
}

export -f get_rep_cls_fun

# run get_rep_cls function using parallel
if [ -s S4_02_AddRepresentID_Cluster.tsv.gz ];then
	echo "SKip. File exists." >&2
else
	pigz -dc -p $threads S4_01_Filter_Cluster_fmt.tsv.gz | parallel -j $threads -k --lb get_rep_cls_fun | pigz -p $threads >S4_02_AddRepresentID_Cluster.tsv.gz

	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> get rep id by parallel: non-zero exit." >&2
		rm -f S4_02_AddRepresentID_Cluster.tsv.gz
		exit 1
	fi
	echo -e "Done. Add rep file:\n$(ls -lh S4_02_AddRepresentID_Cluster.tsv.gz)" >&2
fi

echo "[S4_STEP03]>>> Get the final represent Clustered Sequences ..." >&2
progress_check=$PWD/progress_check
if [ -s "${progress_check}/S4_03.finished" ];then
	echo "Skip running." >&2
	cat ${progress_check}/S4_03.finished >&2
else
	echo ">>>>>> seperate represent IDs by sample to S4_03_EACH_sample_RepresentID dir ..." >&2
	pigz -dc -p $threads S4_02_AddRepresentID_Cluster.tsv.gz | cut -f 1,2 |\
	perl -lane '($sample,$id)=(split(/#/,$F[0]))[0,4];$,="\t";print $sample,$id,"$F[1]#$F[0]";'|\
	csvtk split -j $threads -tTH -f 1 -o S4_03_EACH_sample_RepresentID
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> get split id: non-zero exit." >&2
		rm -rf S4_03_EACH_sample_RepresentID
		exit 1
	fi
	echo "Done." >&2
	echo ">>>>>> Get represent sequences for each sample to S4_04_EACH_sample_RepresentSequences dir ..." >&2
	mkdir -p S4_04_EACH_sample_RepresentSequences
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
	export ttt=$(ls S4_03_EACH_sample_RepresentID/*.tsv | wc -l)
	echo "Progress:( Total records: $ttt )
start: $(date)
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2
	# creat tmp files for counting (use tmp files instead of VARs for FIFO-loop)
	rm -f progress_bar_kkk.tmp progress_bar_ccc.tmp
	echo "start" > progress_bar_ccc.tmp
	touch progress_bar_kkk.tmp
	#####################################
	ls S4_03_EACH_sample_RepresentID/*.tsv | while read cls_file
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
		csvtk join -tTH -f "2;1" -k $cls_file <(seqkit fx2tab $cur_fa) -o S4_04_EACH_sample_RepresentSequences/${tpre}_repSeq.tsv.gz
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
		echo "[ERROR] --> get rep seq tfa: non-zero exit."
		exit 1
	fi
	echo -e "\nDone." >&2
	echo ">>>>>> Merge represent sequences.tfa to S4_04_RepresentSequences_tfa.tsv.gz ..." >&2
	cat S4_04_EACH_sample_RepresentSequences/*_repSeq.tsv.gz > S4_04_RepresentSequences_tfa.tsv.gz
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> merge rep seq: non-zero exit." >&2
		rm -f S4_04_RepresentSequences_tfa.tsv.gz
		exit 1
	fi
	echo ">>>>>> Get final represent sequnces.fasta to S4_04_RepresentSequences.fasta ..." >&2
	csvtk cut -j $threads -tTH -f 3,4 S4_04_RepresentSequences_tfa.tsv.gz | seqkit tab2fx -j $threads -o S4_04_RepresentSequences.fasta
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> get rep seq fasta: non-zero exit." >&2
		rm -f S4_04_RepresentSequences.fasta.gz
		exit 1
	fi
	
	# remove intermediate files if $keep is not "T"
	if [ "$keep" != "T" ];then
		echo "cleaning S4_04_EACH_sample_RepresentSequences and S4_03_EACH_sample_RepresentID ... " >&2
		rm -rf S4_04_EACH_sample_RepresentSequences S4_03_EACH_sample_RepresentID
		echo "Done." >&2
	fi
	echo "Already finished S4_04_EACH_sample_RepresentSequences" >${progress_check}/S4_03.finished
	echo -e "Done. results in: \n$(ls -lh S4_04_RepresentSequences* )" >&2 
fi

# exit ;

echo "
###################################################
 Mapping UNANCHORED seq to Clustered represents
###################################################
"

minimap2 -t $threads -c --cs S4_04_RepresentSequences.fasta S5_UNANCHORED_CLS_rep_seq.fa.gz | pigz -p $threads >S5_UNANCHORED_rep_seq_MM2_S4_Repseq_MergeRef.paf.gz
# echo ">>>> Indexing S4_04_RepresentSequences.fasta ..." >&2
# if [ -s S4_04_RepresentSequences.fasta.bwt ];then
# 	echo "Skip. file exists" >&2
# else
# 	bwa index S4_04_RepresentSequences.fasta
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> bwa index: non-zero exit." >&2
# 		exit 1
# 	fi
# 	echo "DONE." >&2
# fi
# echo ">>>> Get ALL UNANCHORED sequences ..." >&2
# mkdir -p S5_UNANCHORED_Seq
# # get UNANCHORED seq fun
# get_UNANCHORED_seq_fun(){
# 	unanchored_file=$1
# 	if [ ! -s $unanchored_file ];then 
# 		echo "[ERROR]-->No file:$unanchored_file" >&2 
# 		exit 1
# 	fi
# 	unanchored_id=$( cut -f 1 $unanchored_file )
# 	cur_pre=$( basename $unanchored_file )
# 	cur_pre=${cur_pre%%_05_*}
# 	cur_fa=$(grep -m 1 "/${cur_pre}.unaligned" $unaln_seq_list)
# 	if [ ! -s "$cur_fa" ];then
# 		echo "[ERROR] --> No file: $cur_fa" >&2
# 		exit 1
# 	fi
# 	seqkit grep -f <(echo "$unanchored_id") $cur_fa -o S5_UNANCHORED_Seq/${cur_pre}_unanchored.fa.gz
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> get unanchored fa: non-zero exit." >&2
# 		rm -f S5_UNANCHORED_Seq/${cur_pre}_unanchored.fa.gz
# 		exit 1
# 	fi
# 	# run bwa mem for each record
# 	if [ -s "S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq.paf.gz" ];then
# 		echo "Skip MEM. File exists: S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq.paf.gz" >&2
# 	else
# 		bwa mem -t 1 S4_04_RepresentSequences.fasta S5_UNANCHORED_Seq/${cur_pre}_unanchored.fa.gz 2>S5_UNANCHORED_Seq/${cur_pre}_MEM.log | paftools sam2paf - | pigz > S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq.paf.gz
# 		if [[ $? -ne 0 ]] ; then 
# 			echo "[ERROR] --> get MEM PAF: non-zero exit." >&2
# 			rm -f S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq.paf.gz
# 			exit 1
# 		fi
# 	fi
# 	# filter out matched to anchored records and get only unanchored 
# 	echo "Filter for ${cur_pre} with coverage:$coverage and identity:$identity" >&2
# 	zcat S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq.paf.gz |\
# 	perl -lane '$F[0]=~s/\:\.//g;print $F[0] if $F[10]/$F[1] >= $ENV{coverage} && $F[9]/$F[10] >= $ENV{identity}' |\
# 	perl -lane 'print if ++$h{$_}==1' >S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq_filtered.id
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> get filtered id : non-zero exit." >&2
# 		rm -f S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq_filtered.id
# 		exit 1
# 	fi
# 	# get filtered seq
# 	if [ -s "S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq_filtered.id" ];then
# 		seqkit grep -v -f S5_UNANCHORED_Seq/${cur_pre}_MEM_to_RepSeq_filtered.id S5_UNANCHORED_Seq/${cur_pre}_unanchored.fa.gz -o S5_UNANCHORED_Seq/${cur_pre}_unanchored_filtered_out_MEM_to_RepSeq.fa.gz
# 	else
# 		cp S5_UNANCHORED_Seq/${cur_pre}_unanchored.fa.gz S5_UNANCHORED_Seq/${cur_pre}_unanchored_filtered_out_MEM_to_RepSeq.fa.gz
# 	fi
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> get filtered seq: non-zero exit." >&2
# 		rm -f S5_UNANCHORED_Seq/${cur_pre}_unanchored_filtered_out_MEM_to_RepSeq.fa.gz
# 		exit 1
# 	fi
# }
# export -f get_UNANCHORED_seq_fun
# if [ -s "S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz" ];then
# 	echo "Skip. File exists: S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz" >&2
# else
# 	cat $UNANCHORED_list | parallel -j $threads -k --lb get_UNANCHORED_seq_fun 
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> parallel get unanchored: non-zero exit." >&2
# 		exit 1
# 	fi
# 	echo ">>>> Merge each sample unanchored to one ..." >&2
# 	cat S5_UNANCHORED_Seq/*_unanchored_filtered_out_MEM_to_RepSeq.fa.gz > S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> merge unanchored fa: non-zero exit." >&2
# 		rm -f S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz
# 		exit 1
# 	fi
# 	echo "DONE. results in: S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz" >&2
# fi

# echo "
# ###################################################
#  Mapping remained UNANCHORED seq to each other
# ###################################################
# "
# if [ -s "S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz.sa" ];then
# 	echo "Skip index.File exists" >&2
# else
# 	bwa index S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz
# fi
# if [ -s "S5_Filtered_UNANCHORED_MEM_to_self.paf.gz" ];then
# 	echo "Skip MEM. File exists" >&2
# else
# 	bwa mem -t $threads S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz S5_ALL_UNANCHORED_filtered_out_MEM_to_RepSeq_Seq.fa.gz | paftools sam2paf - | pigz -p $threads > S5_Filtered_UNANCHORED_MEM_to_self.paf.gz
# 	if [[ $? -ne 0 ]] ; then 
# 		echo "[ERROR] --> mem to self: non-zero exit." >&2
# 		rm -f S5_Filtered_UNANCHORED_MEM_to_self.paf.gz
# 		exit 1
# 	fi
# 	echo "DONE. results in: S5_Filtered_UNANCHORED_MEM_to_self.paf.gz" >&2 
# fi

echo "
###################################################
Filtering cluster results
###################################################
"

echo "
###################################################
 Merge and get final NON-ref-seq
###################################################
"

