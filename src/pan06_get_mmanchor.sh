#!/usr/bin/env bash

if [ $# -ne 4 ]; then echo "usage $0 <pre> <filtered_unaln.fa> <quast.fcoords> <cpu>" >&2; exit 1;fi

# unaln.fa should have seq-name in this format: XXX_scaffold_123-456

pre=$1
unalnfa=$2
fcoords=$3
cpu=$4
if [ ! -s $fcoords ];then echo "[ERROR]---> No file: $fcoords " >&2;exit 1;fi
if [ ! -s $unalnfa ];then echo "[ERROR]---> No file: $unalnfa " >&2;exit 1;fi
# check format
check_num=$(zcat $unalnfa | head -1 | grep -P "(^.*_scaffold_\d+)_(\d+)-(\d+)"| wc -l )
if [ $check_num -ne 1 ];then 
echo "[ERROR] --> Wrong input fasta format. Please check.
Input fasta should have name in this format:
	[scaffold_id]_[unaln-start]-[unaln-end]
scaffold_id should end with 'scaffold_\d+', eg:
		5A10_scaftig0000000001_scaffold_0_1-286
" >&2;
exit 1;
fi

# get fcoords.bed from fcoords (1-based)
echo "get fcoords.bed1 from $fcoords ..." >&2 &&\
cat $fcoords | perl -lane '$str="+";$s=$F[3];$e=$F[4];if($s>$e){$s=$F[4];$e=$F[3];$str="-";} $,="\t";print $F[12],$s,$e,"$F[11]#$F[0]#$F[1]",$F[9],$str;' >${pre}_fcoords.bed1 &&\
wc -l ${pre}_fcoords.bed1 >&2 &&\
# get unaln_bed from unaln.fa (0-based)
echo "get unaln_bed0 from $unalnfa ..." >&2 &&\
seqkit fx2tab -n $unalnfa | perl -lne 'if(/(^.*_scaffold_\d+)_(\d+)-(\d+)/){$chr=$1;$start=$2-1;$end=$3;print "$chr\t$start\t$end\tunaln\t0\t\+";}' > ${pre}_unaln.bed0 &&\
wc -l ${pre}_unaln.bed0 >&2

if [ -s ${pre}_mmanchor_fmt.bed ];then
	echo "Skip running mmanchor. Output file exist: ${pre}_mmanchor_fmt.bed " >&2
else
	echo "bedtools sort ..." >&2
	if [ ! -s ${pre}_all_sort.bed ];then
		cat ${pre}_fcoords.bed1 ${pre}_unaln.bed0 | bedtools sort >${pre}_all_sort.bed 
		if [[ $? -ne 0 ]] ; then echo "[ERROR] -->get all sort bed: non-zero exit." >&2;rm -f ${pre}_all_sort.bed; exit 1;fi
	else
		echo "Loaded existing sorted bed. Skip sorting." >&2;
	fi
export ttt=$(cat ${pre}_all_sort.bed | grep -c "unaln" );	
	echo "
----------------------------------------------------------------------
Get anchored position by minimap2 alignment coords (filtered in quast).
----------------------------------------------------------------------
Input-minimap-coords: ${pre}_fcoords.bed1
Input-unaln-bed     : ${pre}_unaln.bed0
Sorted-all-bed      : ${pre}_all_sort.bed
Output-mmanchor-bed : ${pre}_mmanchor_fmt.bed
----------------------------------------------------------------------
Output format:
|  1   |   2   |  3  |  4   |   5   |   6    |          7           |         8            |
|------|-------|-----|------|-------|--------|----------------------|----------------------|
| Q_ID | Start | End | Type | Score | Strand | L:chr:posi:str:Qinfo | R:chr:posi:str:Qinfo |

----------------------------------------------------------------------
Progress:( Total records: $ttt )

0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2

echo "start" > progress_bar_ccc.tmp;
rm -f progress_bar_kkk.tmp;
touch progress_bar_kkk.tmp;
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
	# running with multi-threads
	cat ${pre}_all_sort.bed | grep "unaln" | while read ll; 
	do 
	read -u7
	{
		a1l="";b1l="";lls=0;lle=0;a1l_str="";b1l_str="";a1l_posi="";b1l_posi="";
		export lls=$(echo $ll | perl -lane 'print $F[1]');
		export lle=$(echo $ll | perl -lane 'print $F[2]');
		id=$(echo $ll | perl -lane 'print $F[0]'); 
		a1l=$(grep -A1 "$ll" ${pre}_all_sort.bed|tail -1|grep "$id" | grep -v "unaln" );
		b1l=$(grep -B1 "$ll" ${pre}_all_sort.bed|head -1 | grep "$id" | grep -v "unaln" ); 
		if [ -z "$a1l" -a -n "$b1l" ];then 
			# Only left
			type="only-left";
			a1l_posi="NA";
			b1l_posi=$(echo $b1l | perl -lane '($rchr,$rs,$re)=split(/#/,$F[3]);if($F[-1] eq "+"){$posi=$re;}elsif($F[-1] eq "-"){$posi=$rs;} printf "%s:%s:%s:#%s",$rchr,$posi,$F[-1],join(":",@F);')
		elif [ -z "$b1l" -a -n "$a1l" ];then 
			# Only right
			type="only-right";
			a1l_posi=$( echo $a1l | perl -lane '($rchr,$rs,$re)=split(/#/,$F[3]);if($F[-1] eq "+"){$posi=$rs;}elsif($F[-1] eq "-"){$posi=$re;} printf "%s:%s:%s:#%s",$rchr,$posi,$F[-1],join(":",@F);')
			b1l_posi="NA";
		elif [ -n "$a1l" -a -n "$b1l" ];then 
			# Both-end
			type="BothEnd";
			b1l_posi=$(echo $b1l | perl -lane '($rchr,$rs,$re)=split(/#/,$F[3]);if($F[-1] eq "+"){$posi=$re;}elsif($F[-1] eq "-"){$posi=$rs;} printf "%s:%s:%s:#%s",$rchr,$posi,$F[-1],join(":",@F);')
			a1l_posi=$( echo $a1l | perl -lane '($rchr,$rs,$re)=split(/#/,$F[3]);if($F[-1] eq "+"){$posi=$rs;}elsif($F[-1] eq "-"){$posi=$re;} printf "%s:%s:%s:#%s",$rchr,$posi,$F[-1],join(":",@F);')
		elif [ -z "$a1l" -a -z "$b1l" ];then 
			type="Unanchored";
			a1l_posi="NA";
			b1l_posi="NA";
		fi
		lll=$( echo -e "$ll\t$type\tLeft:$b1l_posi\tRight:$a1l_posi" | perl -lane '
		$F[5]=".";
		$F[1]++;
		$len=abs($F[2]-$F[1])+1;
		if(/BothEnd/){
			($lchr,$lstart,$lstr)=(split(/:/,$F[7]))[1,2,3];
			($rchr,$rstart,$rstr)=(split(/:/,$F[8]))[1,2,3];
			if($lchr eq $rchr){ 
				$F[6]=$F[6]."-SameChr"; 
				$ref_len=abs($rstart-$lstart);
				$F[4]=sprintf("%.2f-%s/%s",$ref_len/$len,$ref_len,$len);
				if($lstr eq $rstr){ 
					$F[5]=$lstr; 
					$type="SameStr"; 
					$type2="DiffOrd"; 
					$type2="SameOrd" if $lstr eq "+" and $rstart-$lstart >= 0; 
					$type2="SameOrd" if $lstr eq "-" and $rstart-$lstart < 0; 
					$type=$type."-$type2"; 
				}else{ 
					$type="DiffStr"; 
				} 
				$F[6]=$F[6]."-$type";
			}else{ 
				$F[6]=$F[6]."-DiffChr"; 
			}
		}

		$,="\t";
		print @F[0,1,2,6,4,5,7,8]; ')
		echo "$lll";
		ccc=$(perl -e '$a=`cat progress_bar_ccc.tmp | wc -l`;print int($a*100/$ENV{ttt});')
		kkk=$(cat progress_bar_kkk.tmp | wc -l)
		echo "ccc" >> progress_bar_ccc.tmp
		if [[ $ccc -eq $kkk ]];then
			echo -e "kkk1\nkkk2" >> progress_bar_kkk.tmp
			echo -n "*" >&2;
		fi
		echo >&7 # put the result back to fd
	} &
	done | grep -v "^#" >${pre}_mmanchor_fmt.bed
	wait;
	exec 7>&-    # close fd write
	exec 7<&-    # close fd read
fi
echo -e "\nFinished." >&2
#stats
wc -l  ${pre}_mmanchor_fmt.bed >&2
rm -f progress_bar_ccc.tmp progress_bar_kkk.tmp