#!/usr/bin/env bash
# need third-party programs:
# samtools; csvtk; bedtools; bedops;
if [ $# -ne 5 ]; then echo "usage $0 <pre> <popins:location.txt> <popins:non_ref_new.bam> <mapQ> <cpu>" >&2; exit 1;fi
pre=$1
input=$2
bam=$3
export minq=$4
cpu=$5

if [ -s $input -a -s $bam ];then
	echo "Loaded input files." >&2
else
	echo "[ERROR] --> Can't find $input or $bam !" >&2;
	exit 1;
fi	
# mkdir -p $pre;
# cd $pre;
export r_rg="";
export rstr="";
export qid="";
export qstr="";
export nrp="";
export ac_score="";
export ttt=$(cat $input | wc -l) # progress bar
echo "
Running for $pre
---------------------------------------------------
INPUT: $input
BAM:   $bam
Min MapQ: $minq

OUTPUT:   ${pre}_pploc_qrg_rd_fmt.tsv
------------------------------------------------------
Progress: ( Total records: $ttt )
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2

echo "start" > progress_bar_ccc.tmp;
rm -f progress_bar_kkk.tmp;
touch progress_bar_kkk.tmp;

if [ -s ${pre}_pploc_qrg_rd_fmt.tsv ];then
	echo "Skip running, file exists." >&2
else

	########################################
	# creat pipe for multi-thread running: #
	########################################
	[ -e /tmp/fd1 ] || mkfifo /tmp/fd1	# creat a pipe 
	exec 7<>/tmp/fd1					# creat a File Descriptor (fd) ( read and write )
	rm -rf /tmp/fd1						# rm raw pipe file
	# set thread number to $cpu
	for ((i=1;i<=$cpu;i++))
	do
		echo >&7
	done

	cat $input | while read r_rg rstr qid qstr nrp ac_score
	do
	read -u7
	{
		rsam="";rbed="";read_id="";q_rg="";qbed="";q_read_id="";q_f_rg="";
		#rm -f cur_* ;
		if [[ -z "$r_rg" ]];then
			echo "[ERROR] --> No r_rg" >&2
			echo "[ERROR] --> No r_rg"
			exit 1;
		fi
		if [[ "$r_rg" == "OTHER" ]];then
			echo -e "$r_rg\t$rstr\t$qid\t$qstr\t$nrp\t$ac_score\tUnanchored"
		else
			rsam=$( samtools view $bam $r_rg | grep $qid )
			rbed=$( echo "$rsam" | sam2bed | cut -f 1-6 )
			cur_r_rg_bed=$(echo "$r_rg" | perl -pe 's/[:-]/\t/g')
			read_id=$( echo "$rbed" | intersectBed -f 0.7 -a stdin -b <(echo "$cur_r_rg_bed") | perl -lane 'print $F[3] if $F[4] >= $ENV{minq} and $F[5] eq $ENV{rstr}' )
			q_rg=$( echo "$rsam" | sort -k8,8n | perl -lane '$min=$F[7] if $.==1;$max=$F[7];$id=$F[6];END{print "$id:$min-$max";}')
			qbed=$(samtools view $bam $q_rg | grep -f <( echo "$read_id" ) | sam2bed | cut -f 1-6 | perl -lane 'print if $F[4] >= $ENV{minq} and $F[5] eq $ENV{qstr}')
			q_read_id=$( echo "$qbed" | cut -f 4 | perl -pe 's/\s+/,/g'|perl -pe 's/,$//g' )
			if [ -n "$q_read_id" ];then
				q_f_rg=$( echo "$qbed" | perl -lane '$min=$F[1] if $.==1;$max=$F[2];$id=$F[0];END{print "$id:$min-$max";}' )
				#echo "$qbed" | cut -f 4 >cur_q_id.txt
				r_f_rg=$( echo "$rbed" | grep -f <( echo "$qbed" | cut -f 4 ) | perl -lane '$min=$F[1] if $.==1;$max=$F[2];$id=$F[0];END{print "$id:$min-$max";}' )
				echo -e "$r_f_rg\t$rstr\t$q_f_rg\t$qstr\t$nrp\t$ac_score\t$q_read_id" 
			else
				echo -e "$r_rg\t$rstr\t$qid\t$qstr\t$nrp\t$ac_score\tNo_PASS_Records" 
			fi
		fi
	    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;exit 1;fi
		#-------------progress bar------------------
		ccc=$(perl -e '$a=`cat progress_bar_ccc.tmp | wc -l`;print int($a*100/$ENV{ttt});')
		kkk=$(cat progress_bar_kkk.tmp | wc -l)
		echo "ccc" >> progress_bar_ccc.tmp
		if [[ $ccc -eq $kkk ]];then
			echo -e "kkk1\nkkk2" >> progress_bar_kkk.tmp
			echo -n "*" >&2;
		fi
		#-------------progress bar end---------------
	    echo >&7
	} &
	done | grep -v "^##" > ${pre}_pploc_qrg_rd.txt
	wait;
	exec 7>&-
	exec 7<&-
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_pploc_qrg_rd.txt;exit 1;fi
fi
echo "Formating raw out ${pre}_pploc_qrg_rd.txt to ${pre}_pploc_qrg_rd_fmt.tsv ..."
cat ${pre}_pploc_qrg_rd.txt | perl -lane '
$rid="";$rs="";$re="";$qid="";$qs="";$qe="";
$F[-1]=~/No_PASS|Unanchored/ && next; 
if($F[0]=~/(.*):(\d+)-(\d+)/){ 
	$rid=$1;$rs=$2;$re=$3;
	$rposi=int(($rs+$re)/2); 
}else{
	next;
}
if($F[2]=~/(.*scaffold_\d+)_(\d+)-(\d+).*:(\d+)-(\d+)/){ 
	$qid=$1;$qrs=$2;$qre=$3;$qcs=$4;$qce=$5;
	$qs=$qrs+$qcs;$qe=$qre+$qcs; 
	$qlen_half=int((abs($qre-$qrs)+1)/2); 
	$qposi=int(($qce+$qcs)/2); 
	if($qposi <= $qlen_half){ $type="LEFT"; }else{ $type="RIGHT"; } 
	$,="\t"; 
	print $qid,$qrs,$qre,"PopIns-RD","$qlen_half:$F[4]",$F[3],"$type:$rid:$rposi:$F[1]",@F[0,2,-1]; 
}' | csvtk sort -tTH -j $cpu -k 1 -k 2:n -k 3:n -o ${pre}_pploc_qrg_rd_fmt.tsv
echo -e "\nFinished." >&2
# stats
wc -l ${pre}_pploc_qrg_rd_fmt.tsv >&2
rm -f progress_bar_ccc.tmp progress_bar_kkk.tmp