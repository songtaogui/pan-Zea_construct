#!/usr/bin/env bash
# --- parameters -----
rd1=                    # -1  PE_reads_1.fq.gz
rd2=                    # -2  PE_reads_2.fq.gz
scftig_len=500            # -g  cutoff length for final assembly scaftigs.
cpu=4                    # -t  threads
pre=                    # -x  Prefix for the sample, used to generate outputs.
usr_scaftig=            # -f  will ignore -1 -2 and skip idba_ud
ref=                    # -R  should have bwa indexed
minlen=100                # -l  min length for extract unaligned sequence
blastDB=                # -D  blastDB for remove contaminate. Use NT databse.
plant_nt_id=            # -P  list of all plant id in the blastDB 
tNPSL=0.4               # -N  Non-plant-seq-length-prop cutoff for cleaning scaftigs.
coverage=0.8            # -c  coverage cutoff when filtering primary unalign sequence
identity=0.9            # -i  identity cutoff when filtering primary unalign sequence
raw_bam=                # -B  all reads align to ref raw bam, used to get poorly mapped reads.
read_length=150            # -L  reads length
mapQ=10                    # -q  map Quality cutoff for bam based filtering
tAR=3                    # -A  read-pairs support cutoff 
tSR=3                    # -S  split-reads support cutoff
vcf_header=                # -V  vcf_header file (lines start with "##"), for vcf outputs.
srcdir=                    # -z  Dir for each script. Default: $(dirname $0)/src
# --------------------
workdir=$PWD
srcdir=$(cd $(dirname $0)/src; pwd)
cd $workdir;
# ---------------------

usage=" 
---------------------------------------------------------------------------------------- 
         $(basename "$0") -- identify and anchor non-ref-sequences
---------------------------------------------------------------------------------------- 

# parameters ([R]: Required; [O]: Optional)

-h                         show help and exit.

-1 <file>    [O]    File path of PE_reads_1.fq.gz, coupled with '-2'
-2 <file>    [O]    File path of PE_reads_2.fq.gz, coupled with '-1'
-g <int>     [O]    cutoff length for final assembly scaftigs (Default:500).
-t <int>     [O]    number of threads (Default: 4)
-x <string>  [R]    PREFIX for the sample, used to generate outputs.
-f <file>    [O]    File path of scaftig.fa.gz, will ignore -1 -2 and skip 
                           idba_ud assembly step. If there is properate scaftig 
                           file (Named as: Prefix_scaftig\$scftig_len.fa.gz, eg: 
                           Sample1_scaftig500.fa.gz), the program will use that 
                           file instead of -f scaftig file for downstream analysis. 
                           (NOTE: sequence name of the scaftig file should start 
                           with 'PREFIX_', same as -x PREFIX )
-R <file>    [R]    File path of reference.fa, should have bwa indexed.
-l <int>     [O]    min length for extract unaligned sequences (Default: 100)
-D <file>    [R]    Path of blastDB for cleaning. Use NCBI NT databse.
-P <file>    [R]    Path of list of 'within ids' in the blastDB, one record 
                           per line.(eg.: if your species is plant, list all plant 
                           ids in the blastDB. )
-N <0-1>     [O]    Cleaning scaftig cutoff. Filter out scaftigs with Non-
                           withinSP-seq-length-rate > this value. (Default: 0.4)
-c <0-1>     [O]    coverage cutoff for filtering unaln sequence (Default: 0.8)
-i <0-1>     [O]    identity cutoff for filtering unaln sequence (Default: 0.9)
-B <file>    [R]    File path of all-reads-align-to-ref-raw.bam, used to get 
                           poorly mapped reads in PopIns.
-L <int>     [O]    Input WGS reads length (Default: 150)
-q <int>     [O]    map Quality cutoff for bam based filtering (Default: 10)
-A <int>     [O]    read-pairs support cutoff (Default: 3)
-S <int>     [O]    split-reads support cutoff (Default: 3)
-z <path>    [O]    Path of the sub-scripts dir. (Default: $(dirname $0)/src)
---------------------------------------------------------------------------------------- 
 NOTE: These programs should be put in PATH before running:
    [idba_ud];   [quast5];    [seqkit]; 
    [csvtk];     [pigz];      [blastn]; 
    [bwa];       [samtools];  [bedtools]; 
    [bedops: sam2bed]
    [bbtools: clumpify.sh stats.sh]
    [popins (ignore the velvet program)]
---------------------------------------------------------------------------------------- 
                                                                 Songtao Gui 
                                                                 songtaogui@sina.com 
"

################# GETOPT #########################
if [ $# -eq 0 ]; then echo "$usage"; exit 1;fi

while getopts 1:2:g:t:x:f:R:l:D:P:N:c:i:B:L:q:A:S:V:z:h opt; do 
  case $opt in
    h)  echo "$usage"
        exit
        ;;
    1)  rd1=$OPTARG
        if [ "${rd1:0:1}" != "/" ];then 
            echo "[ERROR] --> [-1]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "$rd1" ]; then
            echo "[ERROR] --> [-1]: No file: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi   
        ;;
    2)  rd2=$OPTARG
        if [ "${rd2:0:1}" != "/" ];then 
            echo "[ERROR] --> [-2]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "$rd2" ]; then
            echo "[ERROR] --> [-2]: No file: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi         
        ;;        
    g)  scftig_len=$OPTARG
        expr $scftig_len "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-g]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2; exit 1      
        fi 
        ;;
    t)  cpu=$OPTARG
        expr $cpu "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-t]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1      
        fi 
        ;;
    x)  pre=$OPTARG
        if [ -z "$pre" ];then
          printf "[ERROR] --> [-x]: PREFIX is empty - %s \n" "$OPTARG" >&2
          echo "$usage" >&2
          exit 1 
        fi
        ;;
    f)  usr_scaftig=$OPTARG
        if [ "${usr_scaftig:0:1}" != "/" ];then 
            echo "[ERROR] --> [-f]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "$usr_scaftig" ]; then
            echo "[ERROR] --> [-f]: No file: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi 
        ;;
    R)  ref=$OPTARG
        if [ "${ref:0:1}" != "/" ];then 
            echo "[ERROR] --> [-R]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "$ref" ]; then
            echo "[ERROR] --> [-R]: No file: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "${ref}.bwt" ];then
            echo "[ERROR] --> [-R]: Ref should be BWA indexed: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi
        ;;
    l)  minlen=$OPTARG
        expr $minlen "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-l]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2; exit 1      
        fi 
        ;;
    D)  blastDB=$OPTARG
        if [ "${blastDB:0:1}" != "/" ];then 
            echo "[ERROR] --> [-D]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [[ ! -s "${blastDB}.00.nhd" ]]; then
            echo "[ERROR] --> [-D]: No file: ${blastDB}*.nhd" >&2;
            echo "$usage" >&2; exit 1
        fi 
        ;;
    P)  plant_nt_id=$OPTARG
        if [ "${plant_nt_id:0:1}" != "/" ];then 
            echo "[ERROR] --> [-P]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "$plant_nt_id" ]; then
            echo "[ERROR] --> [-P]: No file: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi 
        ;;
    N)  tNPSL=$OPTARG
        check_float=$(echo $tNPSL | perl -ne 'print "pass" if $_ > 0 and $_ <= 1;')
        if [ "$check_float" != "pass" ];then
            echo "[ERROR] --> [-N]: illegal number (Should be float between 0-1): $OPTARG" >&2
            echo "$usage" >&2; exit 1
        fi
        ;;
    c)  coverage=$OPTARG
        check_float=$(echo $coverage | perl -ne 'print "pass" if $_ > 0 and $_ <= 1;')
        if [ "$check_float" != "pass" ];then
            echo "[ERROR] --> [-c]: illegal number (Should be float between 0-1): $OPTARG" >&2
            echo "$usage" >&2; exit 1
        fi
        ;;
    i)  identity=$OPTARG
        check_float=$(echo $identity | perl -ne 'print "pass" if $_ > 0 and $_ <= 1;')
        if [ "$check_float" != "pass" ];then
            echo "[ERROR] --> [-i]: illegal number (Should be float between 0-1): $OPTARG" >&2
            echo "$usage" >&2; exit 1
        fi
        ;;
    B)  raw_bam=$OPTARG
        if [ "${raw_bam:0:1}" != "/" ];then 
            echo "[ERROR] --> [-B]: File should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        elif [ ! -s "$raw_bam" ]; then
            echo "[ERROR] --> [-B]: No file: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi
        ;;
    L)  read_length=$OPTARG
        expr $read_length "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-L]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1      
        fi 
        ;;
    q)  mapQ=$OPTARG
        expr $mapQ "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-q]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1      
        fi
        ;;
    A)  tAR=$OPTARG
        expr $tAR "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-A]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1      
        fi
        ;;
    S)  tSR=$OPTARG
        expr $tSR "+" 10 &> /dev/null
        if [ $? -ne 0 ];then
            printf "[ERROR] --> [-S]: illegal number -%s\n" "$OPTARG" >&2
            echo "$usage" >&2
            exit 1      
        fi
        ;;
    # V)  vcf_header=$OPTARG
    #     if [ "${vcf_header:0:1}" != "/" ];then 
    #         echo "[ERROR] --> [-V]: File should use ABSOLUTE PATH: $OPTARG" >&2;
    #         echo "$usage" >&2; exit 1
    #     elif [ ! -s "$vcf_header" ]; then
    #         echo "[ERROR] --> [-V]: No file: $OPTARG" >&2;
    #         echo "$usage" >&2; exit 1
    #     fi 
    #     ;;
    z)  srcdir=$OPTARG
        if [ "${srcdir:0:1}" != "/" ];then 
            echo "[ERROR] --> [-z]: Should use ABSOLUTE PATH: $OPTARG" >&2;
            echo "$usage" >&2; exit 1
        fi
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
################ GETOPT END ##########################

# Test for presence of required arguments
if [[ -z "$pre" || -z "$ref" || -z "$blastDB" || -z "$plant_nt_id" || -z "$raw_bam" ]]
then
    echo "[ERROR] --> One of the required parameters is missing:" >&2
    echo "[ -x -R -D -P -B ]" >&2
    echo $usage >&2
    exit 1
fi
# Test for alternate between -1/-2 and -f 
if [[ -z "$usr_scaftig" ]];then
    # no -f, check -1/-2
    if [[ -z "$rd1" || -z "$rd2" ]];then
        echo "[ERROR] --> One of the following parameters should be set:" >&2
        echo "[-1 & -2] or [ -f ]" >&2
        echo $usage >&2; exit 1
    fi
else
    # check scaftig name and pre name 
    usr_scaftig_pre=$(seqkit fx2tab -n $usr_scaftig | head -1 |  cut -d"_" -f 1 )
    if [ "$usr_scaftig_pre" = "$pre" ];then
        echo "Loaded scaftig sequence file: $usr_scaftig" >&2
    else
        echo "[ERROR] --> INCONSISTENT between scaftig name start '$usr_scaftig_pre' and prefix '$pre'." >&2
        exit 1
    fi
fi


######################### Check Pan-src scripts #####################
# workdir=$PWD
# srcdir=$(cd $(dirname $0)/src; pwd)
# cd $workdir;
num_src=1;
for src_script in pan00_IDBA_assembly.sh pan04_popins_pipe.sh pan08_combine_mm_pp.sh pan01_quast_pre-unaln.sh pan05_pmrc_filter.sh pan09_MMPPSR.sh pan02_blastNT_clean.sh pan06_get_mmanchor.sh pan03_reMEM_filter.sh pan07_get_pploc_qrg_rd.bash
do
    if [ ! -s "$srcdir/$src_script" ]; then
        printf "[ERROR] --> Missing script: %s \n" $srcdir/$src_script >&2
        let num_src++
    fi
done

if [ $num_src -eq 1 ];then
    echo "All src scripts were found in '$srcdir'. Proceeding ..." >&2
else
    echo "Check the Default src dir or set '-z' option to the src script dir." >&2
    exit 1;
fi
########################################################################
## --- third-party-programs ---
# idba_ud
# quast5
# bbtools: clumpify.sh stats.sh
# seqkit
# csvtk
# pigz
# blastn
# bwa
# popins: no velvet
# samtools
# bedtools
# bedops: sam2bed
## ----------------------------
################ Check third-party-programs ########################
num_prg=1;
for program in idba_ud quast.py clumpify.sh stats.sh seqkit csvtk pigz blastn bwa popins samtools bedtools sam2bed
do
    if ! which $program >/dev/null 2>&1 ; then
        printf "[ERROR] --> Missing %-15s Please add to PATH\n" $program >&2
        #echo "[ERROR] --> Missing ${program}. Please add to PATH"
        let num_prg++
    fi
done
if [ $num_prg -eq 1 ];then
    echo "All third-party-programs were found. Proceeding ..." >&2
else
    exit 1;
fi

######################## Main program start #########################
mkdir -p $pre
cd $pre
# idba

echo "
###################################################
########## parse reference genome ... #############
###################################################
" >&2

if [ ! -s "${ref}.sizes" ];then
    seqkit fx2tab -j $threads -n -l ${ref} | perl -lane '$,="\t";print @F[0,-1];' > ${ref}.sizes
    if [ $? -ne 0 ];then echo "get reference sizes failed: Non-zero exit" >&2; rm -f ${ref}.sizes; exit 1;fi
fi

vcf_header=$PWD/${pre}_vcf_header.txt
if [ ! -s "$vcf_header" ];then
    #get vcf header by parsing ref
    cat ${ref}.sizes | perl -F"\t" -lane '
        BEGIN{
            $,="\t";
            print "##fileformat=VCFv4.1";
        }
        print "##contig=<ID=$F[0],length=$F[1]>";
        END{
            $hd_end= <<"END_MESSAGE";
##FILTER=<ID=PASS,Description=\"All filters passed\">
##FILTER=<ID=LowQual,Description="One end supported with Read-pair evidence, but no precise posi">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
END_MESSAGE
            print $hd_end;
        }
    ' > $vcf_header
    if [ $? -ne 0 ];then echo "get vcf header failed: Non-zero exit" >&2; rm -f $vcf_header; exit 1;fi
fi


echo "
###################################################
########## run assembly using idba_ud #############
###################################################
" >&2

if [ -s "${pre}_scaftig${scftig_len}.fa.gz" ];then
    echo "SKIP assembly step. File exists: ${pre}_scaftig${scftig_len}.fa.gz" >&2
else
    if [ -s "$usr_scaftig" ];then
        echo "SKIP assembly step. Using user defined scaftig file: $usr_scaftig" >&2;
    else
        # pan00_IDBA_assembly.sh <pre> <rd1.fq.gz> <rd2.fq.gz> <scaftig-cut-off> <cpus>
        echo "CMD>>> bash $srcdir/pan00_IDBA_assembly.sh $pre $rd1 $rd2 $scftig_len $cpu" >&2
        bash $srcdir/pan00_IDBA_assembly.sh $pre $rd1 $rd2 $scftig_len $cpu
        if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan00_IDBA_assembly.sh: non-zero exit." >&2;exit 1;fi
        # OUTPUTS
        # ${pre}_dedupe.fa.gz
        # idba_out/contig.fa.gz and idba_out/scaffold.fa.gz
        # ${pre}_scaftig${scftig_len}.fa.gz
    fi
fi

echo "
###################################################
## align use quast5 and get primary unalign info ##
###################################################
" >&2
if [ -s "${pre}_scaftig${scftig_len}.fa.gz" ];then
    echo "Found scaftig file: ${pre}_scaftig${scftig_len}.fa.gz" >&2
elif [ -s "$usr_scaftig" ];then
    echo "Using user-set scaftig '$usr_scaftig' as downstream scaftig file:" >&2
    echo "cp -f $usr_scaftig ${pre}_scaftig${scftig_len}.fa.gz" >&2
    cp -f $usr_scaftig ${pre}_scaftig${scftig_len}.fa.gz
fi

scaftig=${pre}_scaftig${scftig_len}.fa.gz
# pan01_quast_pre-unaln.sh <pre> <scaftig.fa> <ref.fa> <minlen> <scftig_len> <cpu>
echo "CMD>>> bash $srcdir/pan01_quast_pre-unaln.sh $pre $scaftig $ref $minlen $scftig_len $cpu" >&2
bash $srcdir/pan01_quast_pre-unaln.sh $pre $scaftig $ref $minlen $scftig_len $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan01_quast_pre-unaln.sh: non-zero exit." >&2;exit 1;fi
# OUTPUTS: 
#  ${pre}_quast dir
#  ${pre}.unaligned.${minlen}bp.bed
#  ${pre}.unaligned.${minlen}bp.fa.gz

# blastNT
echo "
###################################################
 run blastNT and filter Best-Hit-Non-Plant records 
###################################################
" >&2
unaln_info=${pre}_quast/contigs_reports/contigs_report_${scaftig%%.fa*}.unaligned.info
# pan02_blastNT_clean.sh <1:prefix> <2:unaln.fa> <3:scftig.fa> <4:unaln_info> <5:cpus> <6:blastDB> <7:plant_nt_id> <8:NPSL-cutoff> <9:minlen> <10:scftig_len>
echo "CMD>>> bash $srcdir/pan02_blastNT_clean.sh $pre ${pre}.unaligned.${minlen}bp.fa.gz $scaftig $unaln_info $cpu $blastDB $plant_nt_id $tNPSL $minlen $scftig_len" >&2
bash $srcdir/pan02_blastNT_clean.sh $pre ${pre}.unaligned.${minlen}bp.fa.gz $scaftig $unaln_info $cpu $blastDB $plant_nt_id $tNPSL $minlen $scftig_len
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan02_blastNT_clean.sh: non-zero exit." >&2;exit 1;fi
# OUTPUTS:
# ${pre}.blastNTout.gz
# ${pre}_clean_all_scftig${scftig_len}.fa.gz
# ${pre}_clean_unaln${minlen}_scftig${scftig_len}.fa.gz
# ${pre}_unrefseq_kept.fa.gz
# ${pre}_unrefseq_rmNPlant.fa.gz
# ${pre}_rmNPlant.*
# FF_*

# reMEM
echo "
###################################################
### run bwa MEM and filter coverage && identity ###
###################################################
" >&2
# pan03_reMEM_filter.sh <pre> <unrefseq> <MemRef> <coverage> <identity> <cpu>
echo "CMD>>> bash $srcdir/pan03_reMEM_filter.sh $pre ${pre}_unrefseq_kept.fa.gz $ref $coverage $identity $cpu" >&2
bash $srcdir/pan03_reMEM_filter.sh $pre ${pre}_unrefseq_kept.fa.gz $ref $coverage $identity $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan03_reMEM_filter.sh: non-zero exit." >&2;exit 1;fi
# OUTPUTS:
# ${pre}_reMEM.paf.gz
# ${pre}_reMEM_filter.id
# ${pre}.unaligned.MEMfiltered.fa.gz

# popins anchor
echo "
###################################################
########### run popins pipeline ###################
###################################################
" >&2
# pan04_popins_pipe.sh <prefix> <raw_bam> <reMEM_filter_scftg> <reference> <read_length> <cpu> <minlen>
echo "CMD>>> bash $srcdir/pan04_popins_pipe.sh $pre $raw_bam ${pre}.unaligned.MEMfiltered.fa.gz $ref $read_length $cpu $minlen" >&2
bash $srcdir/pan04_popins_pipe.sh $pre $raw_bam $PWD/${pre}.unaligned.MEMfiltered.fa.gz $ref $read_length $cpu $minlen
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan04_popins_pipe.sh: non-zero exit." >&2;exit 1;fi
# OUTPUTS:
# pp_${pre}/${pre} dir

# popins anchor
echo "
###################################################
# filter by poorly-aligned-to-ref reads coverage  #
###################################################
" >&2
popins_bam=pp_${pre}/${pre}/non_ref_new.bam
# pan05_pmrc_filter.sh <pre> <unaln.fa> <popins.bam> <mapQ> <cpu>
echo "CMD>>> bash $srcdir/pan05_pmrc_filter.sh $pre ${pre}.unaligned.MEMfiltered.fa.gz $popins_bam $mapQ $cpu" >&2
bash $srcdir/pan05_pmrc_filter.sh $pre ${pre}.unaligned.MEMfiltered.fa.gz $popins_bam $mapQ $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan05_pmrc_filter.sh: non-zero exit.">&2;exit 1;fi
# OUTPUTS:
# ${pre}_cvrg_dp_mapQ$mapQ.tsv
# ${pre}_pmrcvrg_kept.id
# ${pre}.unaligned.pmrcfiltered.fa.gz


# get anchored info minimap
echo "
###################################################
###### get anchored info from minima2 aln #########
###################################################
" >&2
fcoords=${pre}_quast/contigs_reports/minimap_output/${scaftig%%.fa*}.coords.filtered
# pan06_get_mmanchor.sh <pre> <filtered_unaln.fa> <quast.fcoords> <cpu>
echo "CMD>>> bash $srcdir/pan06_get_mmanchor.sh $pre ${pre}.unaligned.pmrcfiltered.fa.gz $fcoords $cpu" >&2
bash $srcdir/pan06_get_mmanchor.sh $pre ${pre}.unaligned.pmrcfiltered.fa.gz $fcoords $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan06_get_mmanchor.sh: non-zero exit." >&2;exit 1;fi
# OUTPUTS
# ${pre}_fcoords.bed1
# ${pre}_unaln.bed0
# ${pre}_fcoords.bed1
# ${pre}_all_sort.bed
# ${pre}_mmanchor_fmt.bed


# get anchored info read pairs
echo "
###################################################
## get poorly-aligned-to-ref read-group info ######
###################################################
" >&2
# pan07_get_pploc_qrg_rd.bash <pre> <popins:location.txt> <popins:non_ref_new.bam> <mapQ> <cpu>
popins_loc=pp_${pre}/clean_nonrefseq_${minlen}_${pre}_locations.txt
echo "CMD>>> bash $srcdir/pan07_get_pploc_qrg_rd.bash $pre $popins_loc $popins_bam $mapQ $cpu" >&2
bash $srcdir/pan07_get_pploc_qrg_rd.bash $pre $popins_loc $popins_bam $mapQ $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan07_get_pploc_qrg_rd.bash: non-zero exit." >&2;exit 1;fi

# OUTPUTS
# ${pre}_pploc_qrg_rd_fmt.tsv

# combine mm and pp
echo "
###################################################
#### combine align-anchor and read-group info #####
###################################################
" >&2
# pan08_combine_mm_pp.sh <pre> <mmanchor_fmt.bed> <pploc_qrg_rd_fmt.tsv> <cpu>
echo "CMD>>> bash $srcdir/pan08_combine_mm_pp.sh $pre ${pre}_mmanchor_fmt.bed ${pre}_pploc_qrg_rd_fmt.tsv $cpu" >&2
bash $srcdir/pan08_combine_mm_pp.sh $pre ${pre}_mmanchor_fmt.bed ${pre}_pploc_qrg_rd_fmt.tsv $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan08_combine_mm_pp.sh: non-zero exit." >&2;exit 1;fi

# OUTPUTS
# ${pre}_anchor_combine_mm_pp.tsv

# get anchored info : final
echo "
###################################################
###### combine mm pp sr and get final vcf #########
###################################################
" >&2
popins_vcf=pp_${pre}/clean_nonrefseq_${minlen}_${pre}_insertions.vcf
# pan09_MMPPSR.sh <1:pre> <2:AR threshold> <3:SR threshold> <4:vcf_header> <5:popins.vcf> <6:anchor_combine_mm_pp.tsv> <7:cpu>
echo "CMD>>> bash $srcdir/pan09_MMPPSR.sh $pre $tAR $tSR $vcf_header $popins_vcf ${pre}_anchor_combine_mm_pp.tsv $cpu" >&2
bash $srcdir/pan09_MMPPSR.sh $pre $tAR $tSR $vcf_header $popins_vcf ${pre}_anchor_combine_mm_pp.tsv $cpu
if [[ $? -ne 0 ]] ; then echo "[ERROR] --> pan09_MMPPSR.sh: non-zero exit." >&2;exit 1;fi

# clean intimate files
echo "
###################################################
######           cleaning ....            #########
###################################################
" >&2

rm -f ${pre}_anchor_mmpp.tsv ${pre}_anchor_ppSR.tsv ${pre}_MMPP.vcf ${pre}_PPSR.vcf ${pre}_mergeDUP.tsv ${pre}_PPSR_kept.id ${pre}_merged.vcf ${pre}_TTfilter_kept.id

if [ $? -ne 0 ];then echo "remove intimate files failed: Non-zero exit"; exit 1;fi

echo "Congratulations! All done!"

# # OUTPUTS
# ${pre}_anchor_mmpp.tsv
# ${pre}_anchor_ppSR.tsv
# ${pre}_MMPP.vcf
# ${pre}_PPSR.vcf
# ${pre}_mergeDUP.tsv
# ${pre}_PPSR_kept.id
# ${pre}_merged.vcf
# ${pre}_TTfilter_kept.id
# # final outputs
# ${pre}_00_FINAL_ALL.vcf
# ${pre}_01_FINAL_BothEnd.vcf
# ${pre}_02_FINAL_SeFw.vcf
# ${pre}_03_FINAL_SeRv.vcf
# ${pre}_04_FINAL_PASS.vcf
# ${pre}_05_FINAL_UNANCHORED.tsv
# ${pre}_06_FINAL_mmpp.tsv
# ${pre}_07_FINAL_ppSR.tsv
# ${pre}_08_FINAL_anchored_id.txt
