#!/usr/bin/env bash
if [ $# -ne 7 ]; then
	echo -e "USAGE: 
	$0 <1:pre> <2:AR threshold> <3:SR threshold> <4:vcf_header> <5:popins.vcf> <6:anchor_combine_mm_pp.tsv> <7:cpu>
	AR/SR threshold:  filter popins records with less than this number of read-pairs or split-reads 
	mapQ:             for Poorly-mapped-reads-coverage
	" >&2;
	exit 1;
fi
export pre=$1
export num_AR=$2
export num_SR=$3
header=$4				# vcf header file for final vcf outputs
ppSR=$5					# popins vcf output
mmpp=$6					# combine minimap && read-pairs output				
cpu=$7					# cpu


#-------------PROGRAM start ----------------
echo "
>>>>>>>>>>>>>>>>>> ALL start >>>>>>>>>>>>>>>>>>
`date`
" >&2
if [ -s $mmpp ];then
	echo "loaded mmpp-file for $pre" >&2
else
	echo "[ERROR] --> Can't find $mmpp" >&2;
	exit 1;
fi

if [ -s $ppSR ];then
  echo "loaded ppSR-file for $pre" >&2
else
  echo "[ERROR] --> Can't find $ppSR" >&2;
  exit 1;
fi

echo "
>>>>>>>>>>>>>>>>>> STEP 1 >>>>>>>>>>>>>>>>>> 
Filter minimap2 && read-pairs anchor-outputs:
---------------------------------------------------
Criteria:
>>> BOTH-END:
  1. minimap-both-end && same-chr && same-strand && same-order && distance < 1Mb 
  2. minimap-both-end && both end Supported by read-pairs (With in 2k)
  
>>> SINGLE-END:
  1. failed the \"BOTH-END\" criteria above, but one end Supported by read-pairs (With in 2k)
  2. minimap-single-end && Supported by read-pairs (With in 2k)

>>> LOWQUAL:
	minimap2 failed the criteria above:
  	LOWQUAL:BOTH-END && SINGLE-END

>>> UNANCHOR:
	unachored-by minimap2
---------------------------------------------------
INPUT-FILE:   $mmpp 
" >&2

if [ ! -s ${pre}_anchor_mmpp.tsv ];then

	cat $mmpp  |  perl -lane '
	$,="\t";$nmlpr=0;$nplmr=0;
	$queryid=$F[0]."_".$F[1]."-".$F[2];
	$qs=0;
	$qlen=$F[2]-$F[1]+1;
	$qe=$F[2]-$F[1]+1;
	if ( /only-left/ ){
		# # New:
		# # only = SE：
		if(/PLWI2k/){
			# only-left PLWI2k --> SE_left
			print $queryid,$qs,$qe,"SE_left#$F[3]",@F[6,7],@F[4..9] if $nmlpr<1 && $_=~/WI2k/;
		}else{
			print $queryid,$qs,$qe,"SE_left_LOWQUAL#$F[3]",@F[6,7],@F[4..9];	
		}
	}elsif(/only-right/){

		## NEW rule: only = SE
		if(/PRWI2k/){
			# only-right && PRWI2k --> SE_right
			print $queryid,$qs,$qe,"SE_right#$F[3]",@F[6,7],@F[4..9] if $nplmr<1 && $_=~/WI2k/;
		}else{
			print $queryid,$qs,$qe,"SE_right_LOWQUAL#$F[3]",@F[6,7],@F[4..9];	
		}
	}elsif( /BothEnd/ ){
		$mll=$F[6];($mlc,$mlp)=(split(/:/,$mll))[1,2];	
		$mrr=$F[7];($mrc,$mrp)=(split(/:/,$mrr))[1,2];
		if( /SameOrd/ && abs($mrp-$mlp) < 1e6 ){
			# mmBothEnd && SameChr && SameStrand && SameOrder && distance < 1Mb --> BE_mmSameStr
			print $queryid,$qs,$qe,"BE_mmSameStr#$F[3]",@F[6,7],@F[4..9];
		}elsif( /PLWI2k/ ){
			if( /PRWI2k/ ){
				# mmBothEnd && both-end-PP-WI2k --> BE_mmWI2k	
				print $queryid,$qs,$qe,"BE_mmWI2k#$F[3]",@F[6,7],@F[4..9];
			}else{
				# mmBothEnd && single-end-PP-WI2k --> SE
				print $queryid,$qs,$qe,"SE_left#$F[3]",@F[6,7],@F[4..9];
			}
		}else{
			if( /PRWI2k/ ){
				# mmBothEnd && single-end-PP-WI2k --> SE
				print $queryid,$qs,$qe,"SE_right#$F[3]",@F[6,7],@F[4..9];
			}else{
				# mmBothEnd && No pp-support --> LOWQUAL
				print $queryid,$qs,$qe,"BE_LOWQUAL#$F[3]",@F[6,7],@F[4..9];
			}
		}
	}else{
		if($_ =~ /Unanchored/ ){
			#UNANCHOR
			print $queryid,$qs,$qe,"UNANCHOR#$F[3]",@F[6,7],@F[4..9];
		}else{
			#LOWQUAL
			print $queryid,$qs,$qe,"REST_LOWQUAL#$F[3]",@F[6,7],@F[4..9];
		}
	} 
	' > ${pre}_anchor_mmpp.tsv

	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> [mmpp] non-zero exit." >&2;rm -f ${pre}_anchor_mmpp.tsv; exit 1;fi

else
	echo "Skip STEP 1. File exists" >&2
fi

echo "
OUTPUT-FILE:  $PWD/${pre}_anchor_mmpp.tsv
>>>>>>>>>>>>>>>>>> STEP 1 >>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>> Finished >>>>>>>>>>>>>>>>>
" >&2

echo "
>>>>>>>>>>>>>>>>>> STEP 2 >>>>>>>>>>>>>>>>>>
Filter popins_split_reads_anchored results:
---------------------------------------------------
Criteria:
>>> BOTH-END
	same-chr && same-strand && diff-orient ( \"[[ ]]\" or \"]] [[\" )

>>> SINGLE-END
	Supported by no less than $num_AR read-pairs or $num_SR split-reads.
---------------------------------------------------
INPUT-FILE:  $ppSR
" >&2

if [ ! -s "${pre}_anchor_ppSR.tsv" ];then
#   C[5    B9_scaftig0000008788_scaffold_8787_4283-6033  r   :  85   [
#    $1         $2                                        $3     $4   $5
#  ^(.*)   ($ENV{pre}_scaftig.*_scaffold_\d+_\d+-\d+)   ([fr]):(\d+)([\[\]].*)
	grep -P "[fr]:" $ppSR | perl -lane '
	$,="\t";
	if($F[4]=~/^(.*)($ENV{pre}_scaftig.*_scaffold_\d+_\d+-\d+)([fr]):(\d+)([\[\]].*)/){
		$qid=$2;
		$a="$1$5";
		$str=$3;
		$posi=$4;
		$str=~s/\.//;
		$str=~tr/fr/\+\-/;
		if($qid=~/.*_scaffold_\d+_(\d+)-(\d+)/){$ss=$1;
			$ee=$2;
			$hflen=int(($ee-$ss+1)/2);
			$F[2]=$F[2].":$a:$F[-1]:$str";
			$h{$qid}++;
			if($posi <= $hflen){
				$left{$qid}=$F[2].",".$left{$qid};
				$left_posi{$qid}=$posi.",".$left_posi{$qid};
			}
			if($posi >= $hflen){
				$right{$qid}=$F[2].",".$right{$qid};
				$right_posi{$qid}=$posi.",".$right_posi{$qid};
			}
		}
	}
	END{
		foreach $k (sort keys %h ){
			$left_posi{$k}="NA" unless $left_posi{$k}; 
			$right_posi{$k}="NA" unless $right_posi{$k}; 
			$left_posi{$k}=~s/,$//;$right_posi{$k}=~s/,$//;
			$left{$k}=~s/,$//;
			$right{$k}=~s/,$//; 
			unless($left_posi{$k}=~/\d/){
				$left_posi{$k}="NA";
				$left{$k}="NA";
			}
			unless($right_posi{$k}=~/\d/){
				$right_posi{$k}="NA";
				$right{$k}="NA";
			} 
			print $k,"$left_posi{$k}","$right_posi{$k}",$left{$k},$right{$k};
		}
	}' |\
	perl -lane '
	@lp=split(/,/,$F[1]);
	@rp=split(/,/,$F[2]);
	@ll=split(/,/,$F[3]);
	@rr=split(/,/,$F[4]);
	for($l=0;$l<=$#ll;$l++){ 
		($ll_c,$ll_p,$ll_o,$ll_ff,$ll_str)=(split(/:/,$ll[$l]))[0,1,3,4,-1]; 
		$ll_o=~s/.*(\[\[|\]\]).*/$1/;
		$ll_ar=0;
		$ll_ar=$1 if $ll_ff=~/AR=(\d+);/;
		$ll_sr=0;
		$ll_sr=$1 if $ll_ff=~/SR=(\d+);/;
		for($r=0;$r<=$#rr;$r++){
			($rr_c,$rr_p,$rr_o,$rr_ff,$rr_str)=(split(/:/,$rr[$r]))[0,1,3,4,-1];
			$rr_o=~s/.*(\[\[|\]\]).*/$1/;
			$rr_ar=0;
			$rr_ar=$1 if $rr_ff=~/AR=(\d+);/;
			$rr_sr=0;
			$rr_sr=$1 if $rr_ff=~/SR=(\d+);/;
			$,="\t";
			print $F[0],$lp[$l],$rp[$r],"BE_ppSR","Left:$ll[$l]","Right:$rr[$r]" if $ll_c eq $rr_c and $ll_o ne $rr_o;
			if( $ll_ar >= $ENV{num_AR} or $rr_ar >= $ENV{num_AR} or $ll_sr >= $ENV{num_SR} or $rr_sr >= $ENV{num_SR} ){
				#print SINGLE $F[0],$lp[$l],$rp[$r],$ll[$l],$rr[$r] if $ll_c eq "NA" or $rr_c eq "NA";
				print $F[0],$lp[$l],$rp[$r],"SE_ppSR","Left:$ll[$l]","Right:$rr[$r]" if $ll_c eq "NA" or $rr_c eq "NA";
			}
		}
	}' > ${pre}_anchor_ppSR.tsv
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_anchor_ppSR.tsv; exit 1;fi

else
	echo "Skip STEP 2. File exists." >&2
fi 		

echo "
OUTPUTS:     $PWD/${pre}_anchor_ppSR.tsv
>>>>>>>>>>>>>>>>>> STEP 2 >>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>> Finished >>>>>>>>>>>>>>>>>
`wc -l ${pre}_anchor*.tsv`
" >&2

echo "
>>>>>>>>>>>>>>>>>> STEP 3 >>>>>>>>>>>>>>>>>>
-----------------------------------------------------
Generate vcf format outputs and combine mmpp and ppSR
records. (remove ppSR records that had one end within
2K of the same mmpp records.)
-----------------------------------------------------
" >&2
mmtsv=$PWD/${pre}_anchor_mmpp.tsv
pptsv=$PWD/${pre}_anchor_ppSR.tsv
if [[ ! -s $header ]];then echo "[ERROR] --> Cannot find header file: $header" >&2;exit 1;fi 
if [ -s "$mmtsv" -a -s "$pptsv" ];then 
	echo -e "Loaded input files for $pre. " >&2;
else 
	echo -e "Cannot find file: \n$mmtsv\n or \n$pptsv " >&2;
	exit 1;
fi
if [ ! -s  ${pre}_merged.vcf ];then
	# clean raw files if exists.
	rm -f ${pre}_mergeDUP.tsv
	rm -f ${pre}_PPSR_kept.id
	rm -f ${pre}_*vcf
	# MMPP
	echo "Generating MMPP vcf file ..." >&2;
	cat $header > ${pre}_MMPP.vcf;
	cat $mmtsv | perl -lane '
	BEGIN{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ENV{pre}";}
	# Skip unanchored records.
	next if $F[3]=~/^UNANCHOR/;
	$,="\t";
	$vREF="N";
	$vALT="<INS>";
	$vSVMETHOD="MMPP";
	$vSVLEN=$F[2]-$F[1];
	$vQUAL=".";
	$vFILTER="PASS";
	($lc,$lp,$lstr)=(split(/:/,$F[4]))[1,2,3];
	($rc,$rp,$rstr)=(split(/:/,$F[5]))[1,2,3];
	$svfilter=$F[3];$svfilter=~s/#.*//g;
	$vFILTER="LowQual" if $svfilter=~/LOWQUAL/;
	if($svfilter=~/BE_/){
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
		if(/SE_left/){
			#left
			$vCHR=$lc;$vCHR2=$lc;
			$vPOS=$lp;$vEND=$lp;
			$type="L_L";		
			$strand=$lstr;
			$vCT="5to3" if $strand eq "+"; #ctg:f
			$vCT="3to5" if $strand eq "-"; #ctg:r
		}
		if(/SE_right/){
			#right
			$vCHR=$rc;$vCHR2=$rc;
			$vPOS=$rp;$vEND=$rp;
			$type="R_R";	
			$strand=$rstr;	
			$vCT="3to5" if $strand eq "+"; #ctg:r
			$vCT="5to3" if $strand eq "-"; #ctg:f	
		}
	}
	$vID="$svfilter#$vCT#$type#$F[0]#$strand";
	$vINFO="SVTYPE=INS;SVMETHOD=$vSVMETHOD;CHR2=$vCHR2;END=$vEND;SVLEN=$vSVLEN;CT=$vCT;";
	print $vCHR,$vPOS,$vID,$vREF,$vALT,$vQUAL,$vFILTER,$vINFO,"GT\t1/1";
	'>>${pre}_MMPP.vcf

	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_MMPP.vcf; exit 1;fi

	echo "Generated: $PWD/${pre}_MMPP.vcf" >&2;

	#PPSR
	echo -e "\nGenerating PPSR vcf file ..." >&2;
	cat $header >${pre}_PPSR.vcf
	cat $pptsv | perl -lane '
	BEGIN{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ENV{pre}";}
	$,="\t";
	$vREF="N";
	$vALT="<INS>";
	$vSVMETHOD="PPSR";
	# ppSR $F[1] $F[2] indicates the ctg-posi of popins_vcf "[[ctg:<posi>:" regions, the most left or right posi was used if "NOANCHOR"; 
	# So if there were NA in $F[1] and $F[2], the svlen was calculated as the raw unaln-seq length, otherwise as ($F[2]-$F[1]). 
	if($F[1] eq "NA" or $F[2] eq "NA"){
		if($F[0]=~/_(\d+)-(\d+)$/){
			$vSVLEN=$2-$1+1;
		}else{
			next;
		}
	}else{
		$F[1]=1 if $F[1] eq "0"; # change raw 0 to 1 for batch calc length
		$vSVLEN=$F[2]-$F[1]+1;
	}

	$vQUAL=".";
	$vFILTER="PASS";
	($lc,$lp,$lstr)=(split(/:/,$F[4]))[1,2,-1];
	($rc,$rp,$rstr)=(split(/:/,$F[5]))[1,2,-1];
	$svfilter=$F[3];$svfilter=~s/#.*//g;
	if($svfilter=~/BE_/){
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
		## SR SE forward VS reverse:
		# [[ & + : forward  N[[ctg:f
		# [[ & - : forward  N[[ctg:r
		# ]] & + : reverse  ctg:f]]N
		# ]] & - : reverse  ctg:r]]N
		if(/Right:NA/ or /Right:NoPASS/){
			#left
			$vCHR=$lc;$vCHR2=$lc;
			$vPOS=$lp;$vEND=$lp;
			$type="L_L";
			$strand="$lstr";
			$vCT="5to3" if $F[4]=~/\[\[/; 
			$vCT="3to5" if $F[4]=~/\]\]/; 
		}
		if(/Left:NA/ or /Left:NoPASS/){
			#right
			$vCHR=$rc;$vCHR2=$rc;
			$vPOS=$rp;$vEND=$rp;
			$type="R_R";	
			$strand="$rstr";	
			$vCT="5to3" if $F[5]=~/\[\[/; 
			$vCT="3to5" if $F[5]=~/\]\]/; 	
		}
	}
	$vID="$svfilter#$vCT#$type#$F[0]#$strand";
	$vINFO="SVTYPE=INS;SVMETHOD=$vSVMETHOD;CHR2=$vCHR2;END=$vEND;SVLEN=$vSVLEN;CT=$vCT;";
	print $vCHR,$vPOS,$vID,$vREF,$vALT,$vQUAL,$vFILTER,$vINFO,"GT\t1/1" ;
	'>> ${pre}_PPSR.vcf
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_PPSR.vcf; exit 1;fi
	echo -e "Generated: $PWD/${pre}_PPSR.vcf\n\n" >&2;


	duplist=$(cat $mmtsv $pptsv | grep -v "UNANCHOR" | cut -f 1 | sort | uniq -d)
	export ttt=$(echo "$duplist" | wc -l)
	echo "Dealing with $ttt duplacates between $mmtsv and $pptsv ..." >&2
	echo "
Progress: ( Total records: $ttt )
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

	# use hash here?
	# export dupid="";
	echo "$duplist" | while read dupid
	do
	read -u7
	{
		# --- progress bar
		if [ $ttt -gt 50 ];then
			let iii=ccc*100/ttt
			let ccc+=1
			if [[ $iii -eq $kkk ]];then
				#bar=$bar"*"
				let kkk+=2
				echo -n "$bar"
			fi
		fi
		# --- done progress bar
		dupSR_bed=$(
		cat ${pre}_PPSR.vcf |\
		grep $dupid |\
		perl -lane '
		if(/CHR2=(.*?);END=(\d+);/){
			$chr2=$1;
			$end=$2;
			($tt,undef,$oo)=(split(/#/,$F[2]))[0,1,2];
			$oo=~tr/LR/\+\-/;
			($o1,$o2)=split(/_/,$oo);
			$,="\t";
			print @F[0,1],$F[1]+1,$tt,"$F[0]_TtT_$F[1]_TtT_$F[2]",$o1,"SR";
			print $chr2,$end,$end+1,$tt,"$F[0]_TtT_$F[1]_TtT_$F[2]",$o2,"SR";
		}' | uniq
		)

		dupMM_bed=$(
		cat ${pre}_MMPP.vcf |\
		grep $dupid |\
		perl -lane '
		if(/CHR2=(.*?);END=(\d+);/){
			$chr2=$1;
			$end=$2;
			($tt,undef,$oo)=(split(/#/,$F[2]))[0,1,2];
			$oo=~tr/LR/\+\-/;
			($o1,$o2)=split(/_/,$oo);
			$,="\t";
				print @F[0,1],$F[1]+1,"$tt","1",$o1,"MM";
				print $chr2,$end,$end+1,"$tt","1",$o2,"MM";
		}' | uniq
		)

		merge_file=$(
		echo -e "$dupSR_bed\n$dupMM_bed" | sort -k1,1 -k2,2n -k3,3n |\
		mergeBed -i stdin -d 2000 -o collapse -c 1,2,3,4,5,6,7 |\
		sed "s/^/$dupid\t/"
		)
		echo "$merge_file" 
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
	done | grep -v "^#" > ${pre}_mergeDUP.tsv 
	wait;
	exec 7>&-
	exec 7<&-
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_mergeDUP.tsv; exit 1;fi

	grep -v "MM" ${pre}_mergeDUP.tsv | cut -f 9 | sed 's/,/\n/g' | sort -u | sed 's/_TtT_/\t/g' > ${pre}_PPSR_kept.id

	echo -e "\nDONE.\nKept $(wc -l ${pre}_PPSR_kept.id) PPSR records." >&2

	echo "Generating merged vcf file ..." >&2

	cp ${pre}_MMPP.vcf ${pre}_merged.vcf
	grep -f ${pre}_PPSR_kept.id ${pre}_PPSR.vcf >> ${pre}_merged.vcf
else
	echo "Skip STEP 3. File exists." >&2
fi

echo "
OUTPUT:" >&2

ls $PWD/${pre}_merged.vcf

echo "
>>>>>>>>>>>>>>>>>> STEP 3 >>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>> Finished >>>>>>>>>>>>>>>>>
" >&2

echo "
>>>>>>>>>>>>>>>>>> STEP 4 >>>>>>>>>>>>>>>>>
----------------------------------------------------
Filtering vcf with previous MEM and pmrc filter info
----------------------------------------------------
" >&2
if [ ! -s ${pre}_00_FINAL_ALL.vcf ];then
	if [ -s "${pre}_pmrcvrg_kept.id" -a -s "${pre}_reMEM_filter.id" ];then
		echo "Total kept records after filter:" >&2
		csvtk grep -j $cpu -v -P ${pre}_reMEM_filter.id ${pre}_pmrcvrg_kept.id -o ${pre}_TTfilter_kept.id &&\
		wc -l ${pre}_TTfilter_kept.id >&2
		#cat ${pre}_pmrcvrg_kept.id | grep -v -f ${pre}_reMEM_filter.id | sort | uniq -u >${pre}_TTfilter_kept.id
		echo "Generating final outputs ..." >&2
		cat ${pre}_merged.vcf | perl -lane '
		BEGIN{
			open(IN,"$ENV{pre}_TTfilter_kept.id");
			while(<IN>){
				chomp;
				$h{$_}=1;
			}
		}
		/^#/ && print && next;
		$id=(split(/#/,$F[2]))[3];
		print if $h{$id};
		' > ${pre}_00_FINAL_ALL.vcf
	else
		echo "[WARNING] --> No filter info: ${pre}_pmrcvrg_kept.id or ${pre}_reMEM_filter.id , skip filtering. " >&2
		cp ${pre}_merged.vcf ${pre}_00_FINAL_ALL.vcf
	fi
else
	echo "Skip STEP4. File exists." >&2
fi

echo "Generating final results ..." >&2
# seperat BE SE PASS LowQuial
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[2]=~/BE/ && print && next;' >${pre}_01_FINAL_BothEnd.vcf
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[2]=~/SE.*5to3/ && print && next;' >${pre}_02_FINAL_SeFw.vcf
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[2]=~/SE.*3to5/ && print && next;' >${pre}_03_FINAL_SeRv.vcf
cat ${pre}_00_FINAL_ALL.vcf | perl -lane '/^#/ && print && next; $F[6]=~/PASS/ && print && next;' >${pre}_04_FINAL_PASS.vcf
cat ${pre}_00_FINAL_ALL.vcf | grep -v "^#" | perl -lane '$id=(split(/#/,$F[2]))[3];print $id;' > ${pre}_08_FINAL_anchored_id.txt
if [ -s ${pre}_TTfilter_kept.id ];then
	csvtk grep -f 1 -tTH -P ${pre}_TTfilter_kept.id $mmtsv -o ${pre}_06_FINAL_mmpp.tsv
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_0*_FINAL_*; exit 1;fi
	csvtk grep -f 1 -tTH -P ${pre}_TTfilter_kept.id $pptsv -o ${pre}_07_FINAL_ppSR.tsv
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_0*_FINAL_*; exit 1;fi
else
	echo "[WARNING] --> No ${pre}_TTfilter_kept.id, skip generating ${pre}_06_FINAL_mmpp.tsv and ${pre}_07_FINAL_ppSR.tsv" >&2
fi
csvtk grep -f 1 -tTH -v -P ${pre}_08_FINAL_anchored_id.txt ${pre}_06_FINAL_mmpp.tsv -o ${pre}_05_FINAL_UNANCHORED.tsv
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> non-zero exit." >&2;rm -f ${pre}_0*_FINAL_*; exit 1;fi
wc -l ${pre}_0*_FINAL_* | grep -v "total" >&2

echo "
>>>>>>>>>>>>>>>>>> STEP 4 >>>>>>>>>>>>>>>>>>
>>>>>>>>>>>>>>>>> Finished >>>>>>>>>>>>>>>>>
The filterd results:
`ls -1 $PWD/${pre}_0*_FINAL_*`
>>>>>>>>>>>>>>>>>> ALL DONE! >>>>>>>>>>>>>>>>>>
`date`
" >&2
