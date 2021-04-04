#!/usr/bin/env bash
if [ $# -ne 5 ]; then echo "usage $0 <pre> <unaln.fa> <popins.bam> <mapQ> <cpu>" >&2; exit 1;fi

pre=$1
infa=$2
pmrbam=$3
mapQ=$4
cpu=$5

echo "
-----------------------------------------------
Filtering if Poorly-mapped-reads-coverage < 50% 
-----------------------------------------------
" >&2

if [ -s "${pre}.unaligned.pmrcfiltered.fa.gz" -a -s "${pre}_pmrcvrg_kept.id" ];then
	echo "SKIP. Output exists: ${pre}.unaligned.pmrcfiltered.fa.gz" >&2
	seqkit stats -j $cpu $infa ${pre}.unaligned.pmrcfiltered.fa.gz >&2
	exit
fi

if [ -s $pmrbam ];then echo "loaded $pmrbam" >&2;else echo "[ERROR] --> No file $pmrbam" >&2;exit 1;fi
echo "calculating Poorly-mapped-reads-coverage with mapQ = $mapQ ..." >&2
if [ ! -s ${pre}_cvrg_dp_mapQ$mapQ.tsv ] ;then
    echo "get depth and coverage ..." >&2
    samtools depth -Q $mapQ $pmrbam |\
	perl -lane '
	$pos{$F[0]}++;
	$depth{$F[0]}+=$F[2];
	END{
		$,="\t";
		foreach(sort keys %pos){
			print $_,$pos{$_},$depth{$_}/$pos{$_};
		}
	}' | grep "^$pre" |\
	perl -lane '
	if($F[0]=~/scaffold_\d+_(\d+)-(\d+):/){
		$s=$1;
		$e=$2;
	}
	$len=$e-$s+1;
	$,="\t";
	print $F[0],$len,$F[1],$F[1]/$len,$F[2];
	' > ${pre}_cvrg_dp_mapQ$mapQ.tsv
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pmlcvrg: non-zero exit." >&2;exit $?;fi   
else
	echo "Skip. file exists." >&2
fi

echo "Kept reads coverage >= 50% ..." >&2
cat ${pre}_cvrg_dp_mapQ$mapQ.tsv | perl -lane '$F[0]=~s/\:\.//g; print $F[0] if $F[3] >= 0.5' > ${pre}_pmrcvrg_kept.id
wc -l ${pre}_pmrcvrg_kept.id >&2

seqkit replace -j $cpu -p ":. \d+" -r "" $infa | seqkit grep -j $cpu -n -f ${pre}_pmrcvrg_kept.id -o ${pre}.unaligned.pmrcfiltered.fa.gz
# stats
seqkit stats -j $cpu $infa ${pre}.unaligned.pmrcfiltered.fa.gz >&2
