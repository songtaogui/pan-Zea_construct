#!/usr/bin/env bash
if [ ! $# -eq 10 ];then echo "Usage: $0 <1:prefix> <2:unaln.fa> <3:scftig.fa> <4:unaln_info> <5:cpus> <6:blastDB> <7:plant_nt_id> <8:NPSL-cutoff> <9:minlen> <10:scftig_len> " >&2;exit;fi
echo "[`date`] ----------> ALL START" >&2
pre=$1
IN=$2
scftigfa=$3
unalninfo=$4
CPU=$5
DB=$6
export plntid=$7
export tNPSL=$8
tNPSLp=$(echo $tNPSL | perl -ne 'print $_*100')
minlen=$9
scftig_len=$10


if [[ -s "${pre}.blastNTout.gz" && -s "${pre}_clean_all_scftig${scftig_len}.fa.gz" && -s "${pre}_unrefseq_kept.fa.gz" ]];then
    echo "Skip BlastNT filter, files exist.
    $(ls -1 ${pre}.blastNTout.gz ${pre}_clean_all_scftig${scftig_len}.fa.gz ${pre}_unrefseq_kept.fa.gz)" >&2
    exit ;
fi

if [ -s "$plntid" ];then echo "loaded plant nt id files" >&2; else echo "[ERROR] --> No Plant_nt.out file:$plntid" >&2; exit 1; fi  
if [ -s "${DB}.00.nhd" ];then echo "loaded ncbi nt DB files" >&2; else echo "[ERROR] --> No NT Database: ${DB}" >&2; exit 1; fi 

echo -e "Using $CPU cpus \n# remove Blast NT best-hit-not-plant sequences, \n# and output rest to \"$PWD/${pre}_rmNPlant.fa.gz\" ..." >&2
if [ -s $IN ];then

    echo "[`date +"%m-%d %H:%M"`] blastn ..." >&2
    if [ ! -s ${pre}.blastNTout.gz ];then
        blastn -query <(pigz -p $CPU -d -c $IN) -task megablast -db $DB -outfmt 6 -out ${pre}.blastNTout -max_target_seqs 10 -evalue 1e-5 -perc_identity 0.5 -num_threads $CPU
        if [[ $? -ne 0 ]] ; then
            echo "[ERROR] --> blastn returned non-zero exit status - aborting" >&2
            exit 1;
        fi
        pigz -p $CPU ${pre}.blastNTout
    else
        echo "blast output exist.Skip." >&2
    fi    
    echo "[`date +"%m-%d %H:%M"`] get best Plant hits ids" >&2

    if [ ! -s ${pre}.rmNPlant.ID ];then
        zcat ${pre}.blastNTout.gz | perl -lane '$,="\t";print @F[0,1] if ++$a{$F[0]}==1;' | perl -lane 'BEGIN{open(PDB,"$ENV{plntid}");while(<PDB>){chomp;$pdb{$_}=1;}} $,="\t";print @F unless $pdb{$F[1]}' >${pre}.rmNPlant.ID
        if [[ $? -ne 0 ]] ; then
            echo "[ERROR] --> get rmNPlant.ID returned non-zero exit status - aborting" >&2
            exit 1
        fi
    else
        echo "${pre}.rmNPlant.ID exist.Skip." >&2
    fi
    # check if rmNPlant.ID-only has 0 records
    num_rmNPlant=0; [ -s ${pre}.rmNPlant.ID-only ] && num_rnNPlant=$(cat ${pre}.rmNPlant.ID-only | wc -l)
    if [ $num_rmNPlant -gt 0 ];then
        echo "remove best hits not plant seqs: ${pre}.rmNPlant.ID-only" >&2
        cut -f 1 ${pre}.rmNPlant.ID > ${pre}.rmNPlant.ID-only
        echo "get rmNPlant fa files: ${pre}_unrefseq_rmNPlant.fa.gz" >&2
        seqkit grep -j $CPU -f ${pre}.rmNPlant.ID-only $IN -o ${pre}_unrefseq_rmNPlant.fa.gz
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> seqkit grep returned non-zero exit status - aborting" >&2;exit 1;fi
        echo "get kept fa files:  ${pre}_unrefseq_kept.fa.gz" >&2
        seqkit grep -j $CPU -v -f ${pre}.rmNPlant.ID-only $IN -o  ${pre}_unrefseq_kept.fa.gz
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> seqkit grep returned non-zero exit status - aborting" >&2;exit 1;fi
    else
        echo "No non Plant records were found. copy $IN to ${pre}_kept.fa.gz" >&2
        cp -f $IN ${pre}_unrefseq_kept.fa.gz
    fi
else
    echo "[ERROR] --> no $IN" >&2;exit;
fi

npid=${pre}.rmNPlant.ID-only
unalnfa=$IN

if [ $num_rmNPlant -gt 0 ];then
    if [[ -s "$npid" && -s "$unalninfo" && -s "$scftigfa" && -s "$unalnfa" ]];
    then
        echo "Loaded all require files." >&2
    else
        echo "[ERROR] --> Some input files does not exist, please check: 
        $npid
        $unalninfo
        $scftigfa
        $unalnfa
        " >&2
        exit 1
    fi
    echo "Get NoPlantSeqLen for $pre to FF_${pre}_NoPlantSeqLen.tsv ..." >&2
    if [ -s FF_${pre}_NPSL.tsv ];then
        echo "Skipping. File exist." >&2
    else
        cat $npid | perl -F"_" -lane '
            if($#F==4){
                $F[-1]=~s/:.*$//;
                $,="_";
                ($s,$e)=split(/-/,$F[-1]);
                $len=$e-$s+1;
                $id=join("_",@F[0..$#F-1]);
                print "$id\t$len";
            }' | perl -lane '
            BEGIN{
                print "Contig\tNoPlantSeqLen";
            }
            $h{$F[0]}+=$F[1];
            END{
                foreach $k (sort keys %h){ 
                    printf "%s\t%s\n",$k,$h{$k};
                }
            } ' >FF_${pre}_NPSL.tsv
        if [[ $? -ne 0 ]] ; then 
            echo "[ERROR] --> NoPlantSeqLen returned non-zero exit status - aborting" >&2;
            exit 1;
        fi
    fi

    echo "Join unalninfo and NPSL for $pre to FF_${pre}_unaln_NPSL.tsv ..." >&2
    if [ -s FF_${pre}_unaln_NPSL.tsv ];then
        echo "Skipping. File exist." >&2
    else
        csvtk join -j $CPU -tT FF_${pre}_NPSL.tsv $unalninfo > FF_${pre}_unaln_NPSL.tsv
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> csvtk join returned non-zero exit status - aborting" >&2;exit 1;fi
    fi

    echo "Get NPSL prop > $tNPSL ctgID for $pre to FF_${pre}_NPSLgt${tNPSLp}_id.txt ... " >&2
    if [ -s FF_${pre}_NPSLgt${tNPSLp}_id.txt ];then
        echo "Skipping. File exist." >&2
    else
        cat FF_${pre}_unaln_NPSL.tsv|sed '1d'|perl -lane 'print $F[0] if $F[1]/$F[2] > $ENV{tNPSL}' > FF_${pre}_NPSLgt${tNPSLp}_id.txt
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> NPSL > $tNPSL returned non-zero exit status - aborting" >&2;exit 1;fi
    fi

    echo "Get keep_unaln${minlen}_id for $pre to FF_${pre}_keep_unaln${minlen}_id.txt ... " >&2
    if [ -s FF_${pre}_keep_unaln${minlen}_id.txt ];then
        echo "Skipping. File exist." >&2
    else
        seqkit fx2tab -j $CPU -n -i $unalnfa | perl -pe 's/(^.*)_.*/$1/'|sort -u | csvtk grep -j $CPU -tTH -v -P FF_${pre}_NPSLgt${tNPSLp}_id.txt >FF_${pre}_keep_unaln${minlen}_id.txt
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> keep_unaln${minlen} returned non-zero exit status - aborting" >&2;exit 1;fi
    fi

    echo "Get clean_unaln.fa for $pre to ${pre}_clean_unaln${minlen}_scftig${scftig_len}.fa.gz ... " >&2
    if [ -s ${pre}_clean_unaln${minlen}_scftig${scftig_len}.fa.gz ];then
        echo "Skipping. File exist." >&2
    else
        seqkit grep -j $CPU -n -f FF_${pre}_keep_unaln${minlen}_id.txt $scftigfa -o ${pre}_clean_unaln${minlen}_scftig${scftig_len}.fa.gz
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get clean_unaln${minlen}_scftig${scftig_len}.fa returned non-zero exit status - aborting" >&2;exit $?;fi
    fi

    echo "Get clean_all_scftig${scftig_len}.fa for $pre to ${pre}_clean_all_scftig${scftig_len}.fa.gz ... " >&2
    if [ -s ${pre}_clean_all_scftig${scftig_len}.fa.gz ];then
        echo "Skipping. File exist." >&2
    else
        seqkit grep -j $CPU -n -f FF_${pre}_NPSLgt${tNPSLp}_id.txt -v $scftigfa -o ${pre}_clean_all_scftig${scftig_len}.fa.gz
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> clean_all_scftig${scftig_len}.fareturned non-zero exit status - aborting" >&2;exit $?;fi
    fi
else
    echo "No non Plant records were found. copy $scftigfa to ${pre}_clean_all_scftig${scftig_len}.fa.gz" >&2
    cp $scftigfa ${pre}_clean_all_scftig${scftig_len}.fa.gz
fi
echo "---------------------------
Finished
---------------------------
" >&2
