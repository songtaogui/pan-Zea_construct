if [ $# -ne 4 ]; then echo "usage $0 <pre> <mmanchor_fmt.bed> <pploc_qrg_rd_fmt.tsv> <cpu>" >&2 ; exit 1;fi
pre=$1
mmin=$2
ppin=$3
cpu=$4
output=${pre}_anchor_combine_mm_pp.tsv

if [ -s $mmin -a -s $ppin ];then
	echo "loaded input files" >&2
else
	echo -e "[ERROR] --> Can't find inputs:\n  $mmin\nor\n  $ppin !" >&2
	exit 1;
fi

export qid=""
export qs=""
export qe=""
export type=""
export score=""
export qstrand=""
export mmleft=""
export mmright=""

export ttt=$(cat $mmin | wc -l)
echo "
Running $0 for $pre
---------------------------------------------------
minimap_anchor: $mmin
popins_anchor : $ppin

OUTPUT:   $output
------------------------------------------------------
Progress: ( Total records: $ttt )
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2 
echo "start" > progress_bar_ccc.tmp;
rm -f progress_bar_kkk.tmp;
touch progress_bar_kkk.tmp;
# rm -f $output
# touch $output
if [ -s $output ];then
	echo "Skip running, file exist: $output" >&2
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

	cat $mmin | while read qid qs qe type score qstrand mmleft mmright 
	do
	read -u7
	{
		qqq=$(echo -e "$qid\t$qs\t$qe")
		pprecord=$( grep "$qqq" $ppin )
		# extract pp left/right bed and RD
		pplbed=$( echo "$pprecord" | perl -lane '$value=(split(/:/,$F[4]))[1];($tt,$c,$posi,$str)=split(/:/,$F[6]); $a=$F[7]; $F[7]=~s/[:-]/\t/g; $,="\t"; printf "%s\tpp_%s:%s:%s:RD=%s#%s:%s\n",$F[7],$tt,$a,$str,$value,$F[8],$F[5] if $F[6] =~ /^LEFT/' )
		pprbed=$( echo "$pprecord" | perl -lane '$value=(split(/:/,$F[4]))[1];($tt,$c,$posi,$str)=split(/:/,$F[6]); $a=$F[7]; $F[7]=~s/[:-]/\t/g; $,="\t"; printf "%s\tpp_%s:%s:%s:RD=%s#%s:%s\n",$F[7],$tt,$a,$str,$value,$F[8],$F[5] if $F[6] =~ /^RIGHT/' )	
		# merge all pp left/right to one-line
		pplrcd=$( echo "$pplbed" | cut -f 4 | perl -pe 's/\n/,/g' | perl -pe 's/,$//;' )
		pprrcd=$( echo "$pprbed" | cut -f 4 | perl -pe 's/\n/,/g' | perl -pe 's/,$//;' )
		# extract mm left/right, combine with pp left/right, then merge
		mmlbed=$(echo "$mmleft" | perl -pe 's/:/\t/g' | perl -lane '$,="\t";$F[0]="mm_".$F[0];print $F[1],$F[2],$F[2]+1,join(":",@F) unless /unaln|NA/')
		mmrbed=$(echo "$mmright" | perl -pe 's/:/\t/g' | perl -lane '$,="\t";$F[0]="mm_".$F[0];print $F[1],$F[2],$F[2]+1,join(":",@F) unless /unaln|NA/')
		merge_2k_left=$(echo -e "$mmlbed\n$pplbed" | sort -k1,1 -k2,2n -k3,3n |bedtools merge -d 2000 -i stdin -o collapse -c 4 2>/dev/null)
		merge_2k_right=$(echo -e "$mmrbed\n$pprbed" | sort -k1,1 -k2,2n -k3,3n |bedtools merge -d 2000 -i stdin -o collapse -c 4 2>/dev/null)
		# count merge_left right for:
		# num_mm_l num_pp_l num_pm_l
		num_pm_l=$(echo "$merge_2k_left" | grep "mm_" | grep "pp_" |wc -l)
		num_pp_l=$(echo "$merge_2k_left" | grep "pp_" |wc -l)
		num_mm_l=$(echo "$merge_2k_left" | grep "mm_" |wc -l)

		# num_mm_r num_pp_r num_pm_r
		num_pm_r=$(echo "$merge_2k_right" | grep "mm_" | grep "pp_" |wc -l)
		num_pp_r=$(echo "$merge_2k_right" | grep "pp_" |wc -l)
		num_mm_r=$(echo "$merge_2k_right" | grep "mm_" |wc -l)

		# get the type
		if [[ "$num_pm_l" > 0 ]];then
			pltype="PLWI2k";
		elif [[ "$num_pp_l" > 0 && "$num_mm_l" > 0 ]];then
			pltype="PLWO2k";
		elif [[ "$num_pp_l" > 0 ]];then
			pltype="PLONLY";
		else
			pltype="PLNORP";
			pplrcd="pp_LEFT:NA";
		fi

		if [[ "$num_pm_r" > 0 ]];then
			prtype="PRWI2k";
		elif [[ "$num_pp_r" > 0 && "$num_mm_r" > 0 ]];then
			prtype="PRWO2k";
		elif [[ "$num_pp_r" > 0 ]];then
			prtype="PRONLY";
		else
			prtype="PRNORP";
			pprrcd="pp_RIGHT:NA"
		fi
		echo -e "$qid\t$qs\t$qe\t${pltype}#${prtype}#${type}\t$score\t$qstrand\t$mmleft\t$mmright\t$pplrcd\t$pprrcd"
	    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2; rm -f $output; exit 1;fi
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
	done | grep -v "^##" >$output.tmp
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2; rm -f $output.tmp; exit 1;fi
	wait;
	exec 7>&-
	exec 7<&- 
	csvtk sort -j $cpu -tTH -k 1 -k 2:n -k 3:n $output.tmp -o $output
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2; rm -f $output; exit 1;fi
fi
echo -e "\nFinished." >&2
wc -l $output  >&2
rm -f progress_bar_ccc.tmp progress_bar_kkk.tmp $output.tmp