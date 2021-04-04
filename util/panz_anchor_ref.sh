#!/usr/bin/env bash

# run anchor pipe for REF
if [ $# -ne 3 ]; then echo "usage $0 <query-genome.fa.gz> <prefix> <cpu>"; exit 1;fi

# get ump from ref seq
infagz=$1
export pre=$2
cpu=$3
# ref should bwa indexed
ref=$WORK/ref/B73_V4/Zea_mays.AGPv4.dna.toplevel.fa
# vcf header
header=$WORK/b73_vcf_header_INS.txt
# set min unmap seq length (bp)
export minlen=100
#####################Check input######################
if [ ! -s $infagz ];then echo "[ERROR]---> No fasta file: $infagz";exit 1; fi

expr $cpu + 0 &>/dev/null
if [[ $? -ne 0 ]] ; then echo "[ERROR]---> cpu number must be integer!";exit $?;fi

if [ ! -s $header ];then echo "[ERROR]---> No vcf header file: $header ";exit 1;fi

expr $minlen + 0 &>/dev/null
if [[ $? -ne 0 ]] ; then echo "[ERROR]---> minlen number must be integer!";exit $?;fi

if [ ! -s ${ref}.bwt ];then echo "[ERROR]---> No ref file or ref where not bwa indexed: $ref";exit 1;fi

echo "
###################################################
############# get scaftigs 500  ###################
###################################################
"
if [ -s ${pre}_scaftig500.fa.gz ];then
	echo "Skip running. Output file exist: ${pre}_scaftig500.fa.gz "
else
	## get scaftigs 500
	# get gaps --> complement --> filter < 500 --> get fasta
	if [ -s ${infagz%%.gz}.splitcomponent.agp -a -s ${infagz%%.gz}.split.fasta -a -s ${infagz%%.gz}.splitobject.agp ];then
		echo "Skip jcvi split fasta, file exists: $(ls ${infagz%%.gz}.*)";
	else
		jcvi_format fasta gaps --mingap=1 --split --cpus=$cpu $infagz
		if [[ $? -ne 0 ]] ; then echo "[ERROR] --> jcvi split fasta: non-zero exit.";exit $?;fi
	fi
	# $infagz.sizes
	# $infagz.splitcomponent.agp
	# $infagz.split.fasta
	# $infagz.splitobject.agp

	## rename split.fasta to scaftigs:
	cat ${infagz%%.gz}.splitcomponent.agp |perl -lane 'print if $F[4] and $F[4] ne "N";'| perl -lane '$,="\t";print $F[5],join("_",@F[0,1,2],$ENV{pre},sprintf("scaftig%07s",$.))' > ${pre}_scaftig_split_idmap.txt &&\
	csvtk join -tTH -j $cpu -k ${pre}_scaftig_split_idmap.txt  <(seqkit fx2tab ${infagz%%.gz}.split.fasta) | cut -f 2,3 |\
	seqkit tab2fx -j $cpu | seqkit seq -j $cpu --min-len 500 -o ${pre}_scaftig500.fa.gz

	########### raw script to get scaftigs500, found some un-splitted seq, why?? #############
	# jcvi_format fasta gaps --mingap=1 --cpus=$cpu $infagz &&\
	# # get ref.size and sort gaps.bed
	# seqkit fx2tab -j $cpu -n -i -l $infagz | perl -lane '$,="\t";print @F[0,-1];' > ${pre}.size &&\
	# sortBed -i ${infagz%%.gz}.gaps.bed -g ${pre}.size > ${pre}.gaps.bed
	# # get scaftig500.bed
	# complementBed -i ${pre}.gaps.bed -g ${pre}.size | perl -lane '$,="\t";$len=$F[2]-$F[1];$id=sprintf("%s_scaftig%07s",$ENV{pre},$.);print @F,$id,$len if $len >= 500;' >${pre}.scaftig500.bed
	# # get scaftig500.fa.gz
	# # seqkit subseq -j $cpu --bed ${pre}.scaftig500.bed $infagz -o - | seqkit fx2tab -j $cpu | sed 's/:. /_/;s/-/_/;' | seqkit tab2fx -j $cpu -o ${pre}_scaftig500.fa.gz
	# seqkit subseq -j $cpu --bed ${pre}.scaftig500.bed $infagz -o - | seqkit replace -j $cpu -p ":. " -r "_" | seqkit replace -j $cpu -p "-" -r "_" -o ${pre}_scaftig500.fa.gz
	########### raw script to get scaftigs500, found some un-splitted seq, why?? #############

	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get scaftigs 500: non-zero exit.";exit $?;fi
	echo "generated scaftig500 fasta: ${pre}_scaftig500.fa.gz"
	# stats
	seqkit stats -j $cpu ${pre}_scaftig500.fa.gz 
fi

echo "
###################################################
############# align to ref use quast5 #############
###################################################
"
# run quast 
if [ ! -s ${pre}_scaftig500.fa.gz ];then echo "[ERROR]---> No scaftig file: ${pre}_scaftig500.fa.gz";exit 1;fi

unalninfo=${pre}_quast/contigs_reports/contigs_report_${pre}_scaftig500.unaligned.info
fcoords=${pre}_quast/contigs_reports/minimap_output/${pre}_scaftig500.coords.filtered

if [ -s "$unalninfo" -a -s "$fcoords" ];then 
	echo "Skip running. Output file exist: ${pre}_quast ";	
else
	echo "Running Quast 5 ..."
	quast.py ${pre}_scaftig500.fa.gz --debug --eukaryote --fragmented --no-snps --no-plots --no-icarus --no-sv -R $ref --threads $cpu -o ${pre}_quast 
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> quast: non-zero exit.";exit $?;fi
    echo "Quast5 finished nomarlly."
fi

echo "
###################################################
#### get ump seq from quast unalign info ##########
###################################################
"
if [ -s "${pre}.unaligned.${minlen}bp.fa.gz" -a -s "${pre}.unaligned.${minlen}bp.bed" ];then
	echo "Skip running. Output file exist: 
	${pre}.unaligned.${minlen}bp.bed
	${pre}.unaligned.${minlen}bp.fa.gz
	";	
else
	# get_unaligned_fasta.sh ${pre}_scaftig500.fa.gz $unalninfo $minlen ${pre}
	# get > ${minlen} bp bed (0 based):
	if [ -s $unalninfo ];then
		## get unalign percentage and unalign fragment length
	 	echo "get unaligned bed > ${minlen} bp ..."
		sed '1d' $unalninfo| perl -lane '@un=split(/,/,$F[-1]);$,="\t";foreach $uu (@un){ ($s,$e)=split(/-/,$uu);$len=$e-$s+1; print $F[0],$s-1,$e,"unaln",$len,"+" if $len >= $ENV{minlen};}' >${pre}.unaligned.${minlen}bp.bed
	    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get unaln bed: non-zero exit.";exit $?;fi
	else
		echo "[ERROR]---> No $unalninfo !";exit;
	fi
	# get seq by bed:
	if [ -s ${pre}_scaftig500.fa.gz ];then
		echo "get unaligned fasta > ${minlen} bp ..."
		#seqkit subseq -j $cpu --bed ${pre}.unaligned.${minlen}bp.bed ${pre}_scaftig500.fa.gz -o - | seqkit fx2tab -j $cpu | perl -F"\t" -lane '$F[0]=~s/:.*//;$,="\t";print @F;' | seqkit tab2fx -j $cpu -o ${pre}.unaligned.${minlen}bp.fa.gz
		seqkit subseq -j $cpu --bed ${pre}.unaligned.${minlen}bp.bed ${pre}_scaftig500.fa.gz -o - | seqkit replace -j $cpu -p ":.*" -r "" -o ${pre}.unaligned.${minlen}bp.fa.gz
	    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get unaln seq: non-zero exit.";exit $?;fi
	    echo "generated unaligned fasta: ${pre}.unaligned.${minlen}bp.fa.gz "
	else
		echo "[ERROR]---> No ${pre}_scaftig500.fa.gz  !";exit;
	fi
fi
# stats
wc -l ${pre}.unaligned.${minlen}bp.bed
seqkit stats -j $cpu ${pre}.unaligned.${minlen}bp.fa.gz
# outputs: 
## ${pre}.unaligned.${minlen}bp.bed
## ${pre}.unaligned.${minlen}bp.fa.gz
echo "
###################################################
## MEM to REF and filter out by cvrg and Pid ######
###################################################
"
if [ -s ${pre}_reMEM.paf.gz ];then
	echo "Skip running MEM. Output file exist: ${pre}_reMEM.paf.gz "
else
	echo "Aligning unaligned seq back to REF with bwa MEM ..."
	# reMEM to REF and filter
	bwaidx=$ref
	bwa mem -t $cpu $bwaidx ${pre}.unaligned.${minlen}bp.fa.gz | paftools sam2paf - | pigz -p $cpu > ${pre}_reMEM.paf.gz
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> mem and paf: non-zero exit.";exit $?;fi
	echo "generated MEM aln paf: ${pre}_reMEM.paf.gz"
fi

echo "filtering raw unalign records out if coverage > 80% and identity >90% ..."
if [ -s "${pre}.unaligned.${minlen}bp.MEMfiltered.bed" -a -s "${pre}.unaligned.${minlen}bp.MEMfiltered.fa.gz" ];then
	echo "Skip filtering. Output file exist."
else
	zcat ${pre}_reMEM.paf.gz |\
	perl -lane '$F[0]=~s/\:\.//g;print $F[0] if $F[10]/$F[1] >= 0.8 && $F[9]/$F[10] >= 0.9' |\
	perl -lane 'print if ++$h{$_}==1' >${pre}_reMEM_filter.id &&\
	wc -l ${pre}_reMEM_filter.id &&\
	seqkit grep -v -n -f ${pre}_reMEM_filter.id ${pre}.unaligned.${minlen}bp.fa.gz -o ${pre}.unaligned.${minlen}bp.MEMfiltered.fa.gz &&\
	cat ${pre}.unaligned.${minlen}bp.bed | perl -lane '$a=$F[1]+1;$id="$F[0]_$a-$F[2]";$,="\t";print $id,@F' | csvtk grep -j $cpu -v -f 1 -tTH -P ${pre}_reMEM_filter.id | cut -f 2- > ${pre}.unaligned.${minlen}bp.MEMfiltered.bed
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> filter reMEM: non-zero exit.";exit $?;fi
	echo "generated filtered unaligned bed and seq: 
		${pre}.unaligned.${minlen}bp.MEMfiltered.bed
		${pre}.unaligned.${minlen}bp.MEMfiltered.fa.gz
		"
fi
# stats
seqkit stats -j $cpu ${pre}.unaligned.${minlen}bp.fa.gz ${pre}.unaligned.${minlen}bp.MEMfiltered.fa.gz
wc -l ${pre}.unaligned.${minlen}bp.bed ${pre}.unaligned.${minlen}bp.MEMfiltered.bed | grep -v "total"

# get anchored info
echo "
###################################################
############# get anchored info ###################
###################################################
"
if [ -s ${pre}_mmanchor_fmt.bed ];then
	echo "Skip running MEM. Output file exist: ${pre}_mmanchor_fmt.bed "
else
	# use multi-threads
	# get_mm_and_fmt.sh
	if [ ! -s ${pre}.unaligned.${minlen}bp.MEMfiltered.bed ];then echo "[ERROR]---> No file: ${pre}.unaligned.${minlen}bp.MEMfiltered.bed ";exit 1;fi
	if [ ! -s $fcoords ];then echo "[ERROR]---> No file: $fcoords ";exit 1;fi
	cat $fcoords | perl -lane '$str="+";$s=$F[3];$e=$F[4];if($s>$e){$s=$F[4];$e=$F[3];$str="-";} $,="\t";print $F[12],$s,$e,"$F[11]#$F[0]#$F[1]",$F[9],$str;' >${pre}_fcoords.bed

	echo " bedtools sort ..."
	if [ ! -s ${pre}_all_sort.bed ];then
		cat ${pre}_fcoords.bed ${pre}.unaligned.${minlen}bp.MEMfiltered.bed | bedtools sort >${pre}_all_sort.bed ;
		if [[ $? -ne 0 ]] ; then echo "[ERROR] -->get all sort bed: non-zero exit.";rm -f ${pre}_all_sort.bed; exit 1;fi
	else
		echo "Loaded existing sorted bed. Skip sorting.";
	fi
	echo "
----------------------------------------------------------------------
Get anchored position by minimap2 alignment coords (filtered in quast).
----------------------------------------------------------------------
Input-minimap-coords: ${pre}_fcoords.bed
Input-unaln-bed     : ${pre}.unaligned.${minlen}bp.MEMfiltered.bed
Sorted-all-bed      : ${pre}_all_sort.bed
Output-mmanchor-bed : ${pre}_mmanchor_fmt.bed
----------------------------------------------------------------------
Output format:
|  1   |   2   |  3  |  4   |   5   |   6    |          7           |         8            |
|------|-------|-----|------|-------|--------|----------------------|----------------------|
| Q_ID | Start | End | Type | Score | Strand | L:chr:posi:str:Qinfo | R:chr:posi:str:Qinfo |

----------------------------------------------------------------------
"
########### progress bar ###########
####### scripts before loop ########
export ttt=$(cat ${pre}_all_sort.bed | grep -c "unaln" )
echo "Progress:( Total records: $ttt )
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|" >&2
# creat tmp files for counting (use tmp files instead of VARs for FIFO-loop)
echo "start" > progress_bar_ccc.tmp;
rm -f progress_bar_kkk.tmp;
touch progress_bar_kkk.tmp;
#####################################

	#echo "running ..."
	# ttt=$(cat ${pre}_all_sort.bed | grep -c "unaln" );
	# ccc=0;
	# kkk=0;
	touch ${pre}_mmanchor_fmt.bed;

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
		####Put these scripts within loop ####
		ccc=$(perl -e '$a=`cat progress_bar_ccc.tmp | wc -l`;print int($a*100/$ENV{ttt});')
		kkk=$(cat progress_bar_kkk.tmp | wc -l)
		echo "ccc" >> progress_bar_ccc.tmp
		if [[ $ccc -eq $kkk ]];then
			echo -e "kkk1\nkkk2" >> progress_bar_kkk.tmp
			echo -n "*" >&2;
		fi
		#####################################
		echo >&7 # put the result back to fd
	} &
	done | grep -v "^#" >${pre}_mmanchor_fmt.bed
	# 如果不加 grep 直接输出，就会少 $thread 条记录，WHY？？
	wait;
	exec 7>&-    # close fd write
	exec 7<&-    # close fd read
fi
#stats
wc -l  ${pre}_mmanchor_fmt.bed

echo "
###################################################
############# get anchored vcf  ###################
###################################################
"

if [ ! -s ${pre}_mmanchor_fmt.bed ];then echo "No file: ${pre}_mmanchor_fmt.bed";exit 1; fi
cat $header >${pre}_00_FINAL_ALL.vcf;
cat ${pre}_mmanchor_fmt.bed | perl -lane '
$id="$F[0]_$F[1]-$F[2]";
$len=$F[2]-$F[1]+1;
if($F[3]=~/BothEnd/){
	$subtype=(split(/-/,$F[3]))[-1];
	$type="BEref_$subtype";
}else{
	$type="SEref_right" if /only-right/;
	$type="SEref_left" if /only-left/;
	$type="UNANCHOR" if /Unanchored/;
}
$,="\t";
print $id,0,$len,"$type#$F[3]",@F[6,7];
' | csvtk grep -j $cpu -v -f 1 -tTH -P ${pre}_reMEM_filter.id >${pre}_06_FINAL_mmpp.tsv
#if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get 06 out: non-zero exit.";rm -f ${pre}_06_FINAL_mmpp.tsv;exit 1;fi

grep "UNANCHOR" ${pre}_06_FINAL_mmpp.tsv > ${pre}_05_FINAL_UNANCHORED.tsv
#if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get 05 out: non-zero exit.";rm -f ${pre}_05_FINAL_UNANCHORED.tsv;exit 1;fi

cat ${pre}_06_FINAL_mmpp.tsv | perl -lane '
BEGIN{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ENV{pre}";}
# Skip unanchored records.
next if $F[3]=~/^UNANCHOR/;
$,="\t";
$vREF="N";
$vALT="<INS>";
$vSVMETHOD="MMref";
$vSVLEN=$F[2]-$F[1];
$vQUAL=".";
$vFILTER="PASS";
($lc,$lp,$lstr)=(split(/:/,$F[4]))[1,2,3];
($rc,$rp,$rstr)=(split(/:/,$F[5]))[1,2,3];
$svfilter=$F[3];$svfilter=~s/#.*//g;
if($svfilter=~/BEref/){
	# both-end	
	if($lc eq $rc ){	
		# same-chr;
		if($lp <= $rp){
			$vCHR=$lc;$vCHR2=$rc;
			$vPOS=$lp;$vEND=$rp;
			$strand="$lstr$rstr";
			$type="L_R";
		}else{
			$vCHR=$rc;$vCHR2=$lc;
			$vPOS=$rp;$vEND=$lp;
			$strand="$rstr$lstr";
			$type="R_L";
		}
		#print $lc,$start,$end,$type,$F[0],$id;
	}else{
		# diff-chr：smaller chr as CHR, the other as CHR2
		# chrs are all numeric:
		if($lc>0 and $rc>0){
			# chrs are all numeric:
			if($lc <= $rc){
				$vCHR=$lc;$vCHR2=$rc;
				$vPOS=$lp;$vEND=$rp;
				$strand="$lstr$rstr";
				$type="L_R";
			}else{
				$vCHR=$rc;$vCHR2=$lc;
				$vPOS=$rp;$vEND=$lp;
				$strand="$rstr$lstr";
				$type="R_L";
			}
		}else{
			# compare chrs as characters
			if($lc le $rc){
				$vCHR=$lc;$vCHR2=$rc;
				$vPOS=$lp;$vEND=$rp;
				$strand="$lstr$rstr";
				$type="L_R";
			}else{
				$vCHR=$rc;$vCHR2=$lc;
				$vPOS=$rp;$vEND=$lp;
				$strand="$rstr$lstr";
				$type="R_L";
			}    		
		}
	}
	$vCT="5to3";
}else{
	# single-end
	## mm SE forward VS reverse
	# left  &  +   :   forward  N[[ctg:f
	# left  &  -   :   reverse  ctg:r]]N
	# right &  +   :   reverse  ctg:f]]N
	# right &  -   :   forward  N[[ctg:r
	if(/SEref_left/){
		#left
		$vCHR=$lc;$vCHR2=$lc;
		$vPOS=$lp;$vEND=$lp;
		$type="L_L";		
		($lc,$lp,$lstr)=(split(/:/,$F[4]))[1,2,3];
		$strand="$lstr";
		$vCT="5to3" if $strand eq "+"; #ctg:f
		$vCT="3to5" if $strand eq "-"; #ctg:r
	}elsif(/SEref_right/){
		#right
		$vCHR=$rc;$vCHR2=$rc;
		$vPOS=$rp;$vEND=$rp;
		$type="R_R";	
		($rc,$rp,$rstr)=(split(/:/,$F[5]))[1,2,3];
		$strand="$rstr";	
		$vCT="3to5" if $strand eq "+"; #ctg:r
		$vCT="5to3" if $strand eq "-"; #ctg:f	
	}
}
$vID="$svfilter#$vCT#$type#$F[0]#$strand";
$vINFO="SVTYPE=INS;SVMETHOD=$vSVMETHOD;CHR2=$vCHR2;END=$vEND;SVLEN=$vSVLEN;CT=$vCT;";
print $vCHR,$vPOS,$vID,$vREF,$vALT,$vQUAL,$vFILTER,$vINFO,"GT\t1/1" unless $svfilter=~/UNANCHOR/;
' >>${pre}_00_FINAL_ALL.vcf
# seperat BE SE PASS LowQuial
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[2]=~/BE/ && print && next;' >${pre}_01_FINAL_BothEnd.vcf
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[2]=~/SE.*5to3/ && print && next;' >${pre}_02_FINAL_SeFw.vcf
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[2]=~/SE.*3to5/ && print && next;' >${pre}_03_FINAL_SeRv.vcf
echo "DONE!
$(wc -l ${pre}_0* | grep -v "total")
"
cp ${pre}.unaligned.${minlen}bp.MEMfiltered.fa.gz ${pre}.unaligned.pmrcfiltered.fa.gz

echo "
###################################################
############## ALL FININSHED !  ###################
###################################################
"
