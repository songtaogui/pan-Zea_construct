if [ $# -ne 5 ]; then echo -e "USAGE: $0 <pre> <rd1.fq.gz> <rd2.fq.gz> <scaftig-cut-off: 500bp> <cpus>\n" >&2; exit 1;fi

sample=$1
rd1=$2
rd2=$3
export cutoff=$4
cpu=$5

if [ -s $rd1 -a -s $rd2 ];then
  echo -e "\n--> Reads loaded.\n" >&2
else
  echo "[ERROR] --> Can't find reads for $sample !" >&2
  exit 1
fi

if [[ ! -d ${sample}_idba ]]
then
  mkdir -p ${sample}_idba
fi
cd ${sample}_idba
echo "
==========================
 DEALING WITH: $sample
==========================
" >&2

start=`date +%s`

# dedupe
if [ ! -f "idba_out/end" ];then
  echo "[`date +"%m-%d %H:%M"`] ----> dedupe START" >&2
  #echo "[`date +"%m-%d %H:%M"`] ----> dedupe START" >>$logf
  if [ -s "${sample}_dedupe.fa.gz" ];then
    pigz -p $cpu -d ${sample}_dedupe.fa.gz
  elif [ ! -s "${sample}_dedupe.fa" ];then
    clumpify.sh -Xmx40g threads=$cpu in1=$rd1 in2=$rd2 dedupe subs=0 out=${sample}_dedupe.fa >&2
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> clumpify: non-zero exit." >&2;exit 1;fi
  else
    echo "Skip dedupe. using pre-exist files instead." >&2
  fi
  # assembly
  echo "[`date +"%m-%d %H:%M"`] ----> assembly START" >&2
  #echo "[`date +"%m-%d %H:%M"`] ----> assembly START" >>$logf
    idba_ud --num_threads $cpu --mink 20 --maxk 100 --step 20 -r ${sample}_dedupe.fa -o idba_out >&2
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> idba_ud: non-zero exit." >&2;exit 1;fi
else
  echo "Skip assembly. Already DONE!" >&2
fi

# clean assembly
echo "[`date +"%m-%d %H:%M"`] ----> Clean START" >&2
#echo "[`date +"%m-%d %H:%M"`] ----> Clean START" >>$logf
cd idba_out
if [ -f "end" ];then
  echo "rm idba intermediate files" >&2
  rm -f align-* contig-* graph-* local-contig-* kmer
  echo "compressing results ..." >&2
  if [ -s "contig.fa" -a ! -s "contig.fa.gz" ];then
    pigz -p $cpu -9 contig.fa
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pigz: non-zero exit." >&2;exit 1;fi
  fi
  if [ -s "scaffold.fa" -a ! -s "scaffold.fa.gz" ];then 
    pigz -p $cpu -9 scaffold.fa
    if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pigz: non-zero exit." >&2;exit 1;fi
  fi
fi
# compress dedupe.fa
cd ..  
echo "compressing deduped reads ..." >&2
if [ -s "${sample}_dedupe.fa" ];then 
  pigz -p $cpu -9 ${sample}_dedupe.fa;
  if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pigz: non-zero exit." >&2;exit 1;fi
fi

# get >= cutoff bp scaftig.fa
export sp=$sample
if [ ! -s "${sample}_scaftig${cutoff}.fa.gz" ];then
  echo "get >= $cutoff bp scaftig.fa for $sp" >&2
  #seqkit seq -m ${cutoff}idba_out/scaffold.fa.gz -o ${sample}_ctg${cutoff}.fa.gz 1>>$logf 2>&1
  seqkit fx2tab -j $cpu idba_out/scaffold.fa.gz |\ 
  perl -F"\t" -lane '$F[1]=~s/n+/\n$F[0]\t/gi;$,="\t";print @F' |\ 
  perl -F"\t" -lane '$scftig=sprintf("scaftig%010d",$.);$F[0]="$ENV{sp}_$scftig"."_$F[0]";$,="\t";print @F if length($F[1]) >= $ENV{cutoff};' |\ 
  seqkit tab2fx -j $cpu -o ${sample}_scaftig${cutoff}.fa.gz
  if [[ $? -ne 0 ]] ; then echo "[ERROR] --> get scaftigs: non-zero exit." >&2;exit 1;fi
fi

end=`date +%s`
runtime=$((end-start))
let runtime_m=runtime/60
let runtime_h=runtime/3600
echo "[`date +"%m-%d %H:%M"`] -----[RUNTIME: $runtime s => $runtime_m m => $runtime_h h ]-----> ALL DONE" >&2
