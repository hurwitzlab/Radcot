#!/bin/bash

#SBATCH -J cntrfge 
#SBATCH -A iPlant-Collabs 
#SBATCH -N 4
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>
# Second author: Scott G. Daniel <scottdaniel@email.arizona.edu>

module load tacc-singularity 
module load launcher

set -u

#
# Set up defaults for inputs, constants
#
IN_DIR=""
QUERY=""
MODE="single"
FASTX=""
FORWARD=""
REVERSE=""
SINGLETONS=""
INDEX="p_compressed+h+v"
OUT_DIR="$PWD/centrifuge-out"
INDEX_DIR="/work/05066/imicrobe/iplantc.org/data/centrifuge-indexes"
MAX_SEQS_PER_FILE=1000000
CENTRIFUGE_IMG="centrifuge-patric.img"
EXCLUDE_TAXIDS=""
SKIP_EXISTING=1
#If you have your own launcher setup on stampede2 just point MY_PARAMRUN at it
#this will override the TACC_LAUNCHER...
PARAMRUN="${MY_PARAMRUN:-$TACC_LAUNCHER_DIR/paramrun}"
MIN_ABUNDANCE=0.01
FORMAT="fasta"

#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    printf "Usage:\n  %s -q DIR_OR_FILE\n\n" "$(basename "$0")"
    printf "Usage:\n  %s -d IN_DIR\n\n" "$(basename "$0")"
    printf "Usage:\n  %s -a FASTX\n\n" "$(basename "$0")"
    printf "Usage:\n  %s -1 FASTX_r1 -2 FASTX_r2 [-s SINGLETONS]\n\n" \
      "$(basename "$0")"
  
    echo "Required arguments:"
    echo " -q DIR_OR_FILE"
    echo ""
    echo "OR"
    echo " -d IN_DIR (single-only)"
    echo ""
    echo "OR"
    echo " -a FASTX (single)"
    echo ""
    echo "OR"
    echo " -1 FASTX_r1 (forward)"
    echo " -2 FASTX_r2 (reverse)"
    echo ""
    echo "Options:"
    echo " -f FORMAT ($FORMAT)"
    echo " -i INDEX ($INDEX)"
    echo " -o OUT_DIR ($OUT_DIR)"
    echo " -s SINGLETONS"
    echo " -k SKIP_EXISTING ($SKIP_EXISTING)"
    echo " -m MIN_ABUNDANCE ($MIN_ABUNDANCE)"
    echo " -x EXCLUDE_TAXIDS"
    echo ""
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && HELP

while getopts :a:d:i:f:m:o:q:r:s:m:x:k:1:2:h OPT; do
    case $OPT in
        a)
            FASTX="$OPTARG"
            ;;
        d)
            IN_DIR="$OPTARG"
            ;;
        h)
            HELP
            ;;
        i)
            INDEX="$OPTARG"
            ;;
        f)
            FORMAT="$OPTARG"
            ;;
        k)
            SKIP_EXISTING=1
            ;;
        m)
            MODE="$OPTARG"
            ;;
        o)
            OUT_DIR="$OPTARG"
            ;;
        q)
            QUERY="$QUERY $OPTARG"
            ;;
        1)
            FORWARD="$OPTARG"
            ;;
        2)
            REVERSE="$OPTARG"
            ;;
        s)
            SINGLETONS="$OPTARG"
            ;;
        m)
            MIN_ABUNDANCE="$OPTARG"
            ;;
        x)
            EXCLUDE_TAXIDS="$OPTARG"
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            exit 1
            ;;
        \?)
            echo "Error: Invalid option: -${OPTARG:-""}"
            exit 1
    esac
done

if [[ ! -e "$CENTRIFUGE_IMG" ]]; then
    echo "Missing CENTRIFUGE_IMG \"$CENTRIFUGE_IMG\""
    exit 1
fi

#
# Verify existence of INDEX_DIR, chosen INDEX
#
if [[ ! -d "$INDEX_DIR" ]]; then
    echo "Cannot find INDEX_DIR \"$INDEX_DIR\""
    exit 1
fi

NUM=$(find "$INDEX_DIR" -name $INDEX.\*.cf | wc -l | awk '{print $1}')

if [[ $NUM -gt 0 ]]; then
    echo "Using INDEX \"$INDEX\""
else
    echo "Cannot find INDEX \"$INDEX\""
    echo "Please choose from the following:"
    find "$INDEX_DIR" -name \*.cf -exec basename {} \; | sed "s/\.[0-9]\.cf//" | sort | uniq | cat -n
    exit 1
fi

#
# Verify existence of various directories, files
#
[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

REPORT_DIR="$OUT_DIR/reports"
[[ ! -d "$REPORT_DIR" ]] && mkdir -p "$REPORT_DIR"

PLOT_DIR="$OUT_DIR/plots"
[[ ! -d "$PLOT_DIR" ]] && mkdir -p "$PLOT_DIR"

if [[ ! -d "$TACC_LAUNCHER_DIR" ]]; then
    echo "Cannot find TACC_LAUNCHER_DIR \"$TACC_LAUNCHER_DIR\""
    exit 1
fi

if [[ ! -f "$PARAMRUN" ]]; then
    echo "Cannot find PARAMRUN \"$PARAM_RUN\""
    exit 1
fi

#
# Create, null-out command file for running Centrifuge
#
CENT_PARAM="$PWD/$$.centrifuge.param"
cat /dev/null > "$CENT_PARAM"

EXCLUDE_ARG=""
[[ -n "$EXCLUDE_TAXIDS" ]] && EXCLUDE_ARG="--exclude-taxids $EXCLUDE_TAXIDS"
RUN_CENTRIFUGE="CENTRIFUGE_INDEXES=$INDEX_DIR singularity run $CENTRIFUGE_IMG $EXCLUDE_ARG"

#
# Set up LAUNCHER env
#
#export LAUNCHER_DIR="$HOME/src/launcher"
#export LAUNCHER_PLUGIN_DIR="$LAUNCHER_DIR/plugins"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_RMI=SLURM
export LAUNCHER_SCHED=interleaved

INPUT_FILES=$(mktemp)

#
# A single FASTX
#
if [[ -n "$FASTX" ]]; then
    BASENAME=$(basename "$FASTX")
    echo "Will process single FASTX \"$BASENAME\""
    echo "$RUN_CENTRIFUGE -f -x $INDEX -U $FASTX -S $REPORT_DIR/$BASENAME.sum --report-file $REPORT_DIR/$BASENAME.tsv" > "$CENT_PARAM"

#
# Paired-end FASTX reads
#
elif [[ -n "$FORWARD" ]] && [[ -n "$REVERSE" ]]; then
    BASENAME=$(basename "$FORWARD")
    echo "Will process FORWARD \"$FORWARD\" REVERSE \"$REVERSE\""

    S=""
    [[ ! -z $SINGLETONS ]] && S="-U $SINGLETONS"

    echo "$RUN_CENTRIFUGE -f -x $INDEX -1 $FORWARD -2 $REVERSE $S -S $REPORT_DIR/$BASENAME.sum --report-file $REPORT_DIR/$BASENAME.tsv" > "$CENT_PARAM"

#
# A directory of single FASTX files
#
elif [[ -n "$IN_DIR" ]] && [[ -d "$IN_DIR" ]]; then
    if [[ $MODE == 'single' ]]; then
        find "$IN_DIR" -type f -size +0c \( -name \*.fa -o -name \*.fasta \) > "$INPUT_FILES"
    else
        echo "Can't yet run IN_DIR with 'paired' mode"
        exit 1
    fi

#
# Either files and/or directories
#
elif [[ -n "$QUERY" ]]; then
    for QRY in $QUERY; do
        if [[ -d "$QRY" ]]; then
            find "$QRY" -type f -not -name .\* >> "$INPUT_FILES"
        elif [[ -f "$QRY" ]]; then
            echo "$QRY" >> "$INPUT_FILES"
        else 
            echo "QUERY ARG \"$QRY\" is neither dir nor file"
        fi
    done

#
# Else "error"
#
else
    echo "Must have -q FILE_OR_DIR/-d IN_DIR/-a FASTX/-f FORWARD & -r REVERSE [-s SINGLETON]"
    exit 1
fi

NUM_INPUT=$(lc "$INPUT_FILES")
if [[ $NUM_INPUT -gt 0 ]]; then
    SPLIT_DIR="$OUT_DIR/split"
    [[ ! -d "$SPLIT_DIR" ]] && mkdir -p "$SPLIT_DIR"
    SPLIT_PARAM="$$.split.param"

    i=0
    while read -r FILE; do
        BASENAME=$(basename "$FILE")
        FILE_SPLIT_DIR="$SPLIT_DIR/$BASENAME"
        NUM_SPLIT_FILES=0
        if [[ -d "$FILE_SPLIT_DIR" ]]; then
            NUM_SPLIT_FILES=$(find "$FILE_SPLIT_DIR" -type f | wc -l | awk '{print $1}')
        fi

        if [[ $NUM_SPLIT_FILES -lt 1 ]]; then
            let i++
            printf "%6d: Split %s\n" $i "$(basename "$FILE")"
            echo "singularity exec $CENTRIFUGE_IMG fxsplit.py -i $FILE -f $FORMAT -o $FILE_SPLIT_DIR -n $MAX_SEQS_PER_FILE" >> "$SPLIT_PARAM"
        fi
    done < "$INPUT_FILES"

    echo "Launching splitter"
    export LAUNCHER_PPN=8
    export LAUNCHER_JOB_FILE="$SPLIT_PARAM"
    "$TACC_LAUNCHER_DIR/paramrun"
    rm "$SPLIT_PARAM"

    SPLIT_FILES=$(mktemp)
    find "$SPLIT_DIR" -type f -size +0c > "$SPLIT_FILES"
    NUM_SPLIT=$(lc "$SPLIT_FILES")
    echo "Splitter done, found NUM_SPLIT \"$NUM_SPLIT\""

    while read -r FILE; do
        BASENAME=$(basename "$FILE")
        SUM_FILE="$REPORT_DIR/$BASENAME.sum"
        TSV_FILE="$REPORT_DIR/$BASENAME.tsv"
  
        if [[ "$SKIP_EXISTING" -gt 0 ]] && [[ -s "$SUM_FILE" ]] && [[ -s "$TSV_FILE" ]]; then
            echo "Skipping $BASENAME - sum/tsv files exist"
        else
            if [[ "$FORMAT" == "fasta" ]]; then
    #            echo "This is a fasta" #debug
                echo "$RUN_CENTRIFUGE -f -x $INDEX -U $FILE -S $REPORT_DIR/$BASENAME.sum --report-file $REPORT_DIR/$BASENAME.tsv" >> "$CENT_PARAM"
            elif [[ "$FORMAT" == "fastq" ]]; then
    #            echo "This is a fastq" #debug
                echo "$RUN_CENTRIFUGE -x $INDEX -U $FILE -S $REPORT_DIR/$BASENAME.sum --report-file $REPORT_DIR/$BASENAME.tsv" >> "$CENT_PARAM"
            else
                echo "File is not fasta or fastq!"
                exit 1
            fi
        fi
    done < "$SPLIT_FILES"

    rm "$SPLIT_FILES"
    #rm -rf "$SPLIT_DIR"
fi

#
# Pass Centrifuge run to LAUNCHER
# Run "interleaved" to ensure this finishes before bubble
#
NUM_CENT_JOBS=$(lc "$CENT_PARAM")
if [[ "$NUM_CENT_JOBS" -gt 0 ]]; then
    echo "Running \"$NUM_CENT_JOBS\" for Centrifuge \"$CENT_PARAM\""
    export LAUNCHER_JOB_FILE="$CENT_PARAM"
    export LAUNCHER_PPN=4
    "$LAUNCHER_DIR/paramrun"
    echo "Finished Centrifuge"
else
    echo "There are no Centrifuge jobs to run!"
    exit 1
fi

rm "$CENT_PARAM"
#
# Collapse the results
#
COLLAPSE_DIR="$OUT_DIR/collapsed"
echo "Collapsing reports"
#echo "DEBUG"
#echo "These are the input files: "$INPUT_FILES""
#echo "This is the report dir: "$REPORT_DIR""
#echo "This is the collapse dir: "$COLLAPSE_DIR""
singularity exec $CENTRIFUGE_IMG collapse.py -l "$INPUT_FILES" -r "$REPORT_DIR" -o "$COLLAPSE_DIR"
echo "Finished collapse"

#rm "$INPUT_FILES"

#
# Create bubble plot
#
echo "Starting bubble"
singularity exec $CENTRIFUGE_IMG centrifuge_bubble.r --dir "$COLLAPSE_DIR" --outdir "$PLOT_DIR" --outfile "bubble" --title "centrifuge"

#BUBBLE_PARAM="$PWD/$$.bubble.param"
#echo "singularity exec $CENTRIFUGE_IMG centrifuge_bubble.r --dir $COLLAPSE_DIR --outdir $PLOT_DIR --outfile bubble --title centrifuge" > "$BUBBLE_PARAM"
#export LAUNCHER_JOB_FILE="$BUBBLE_PARAM"
#"$LAUNCHER_DIR/paramrun"
echo "Finished bubble"

#
# Getting genomes from PATRIC
#
GENOME_DIR="$OUT_DIR/genomes"
[[ ! -d $GENOME_DIR ]] && mkdir -p $GENOME_DIR
echo "Getting genomes and annotations from patricbrc.org"
#-r directory with tsv report files -o output directory for genomes and annotations
singularity exec $CENTRIFUGE_IMG cfuge_to_genome.py -r "$COLLAPSE_DIR" -o $GENOME_DIR -m $MIN_ABUNDANCE

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Ken Youens-Clark <kyclark@email.arizona.edu>"
echo "or Scott Daniel <scottdaniel@email.arizona.edu>"
echo "for version that runs patric-cli along with centrifuge (RADCOT)"
#!/bin/bash

#SBATCH -J bowtie2
#SBATCH -A iPlant-Collabs 
#SBATCH -N 12
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal

# Author: Scott G. Daniel <scottdaniel@email.arizona.edu>

###Uncomment when back on tacc#
#echo "#### Current modules after app.json processing:"
#module list 2>&1
echo "#### LOADING TACC-SINGULARITY ####"
module load tacc-singularity 2>&1
echo "#### LOADING LAUNCHER ####"
module load launcher 2>&1
#echo "#### Current modules after run.sh processing:"
#module list 2>&1
#
# Set up defaults for inputs, constants
#
SING_IMG="bowtie_sam.img"
OUT_DIR="$PWD/bowtie-samtools-out"
#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    singularity exec $SING_IMG patric_bowtie2.py -h
    exit 0
}

function build_opt_string() {
    if [[ -n "$2" ]]; then
        OPTSTRING="$OPTSTRING "$1" "$2""
    fi
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

#set -u

#parse options
while getopts :g:x:1:2:U:f:O:n:l:a:e:L:N5:3:I:X:t:A:h ARG; do
    case $ARG in
        g)
            GENOME_DIR="$OPTARG"
            ;;
        x)
            BT2_IDX="$OPTARG"
            ;;
        1)
            M1="$OPTARG"
            ;;
        2)
            M2="$OPTARG"
            ;;
        U)
            UNPAIRED="$OPTARG"
            ;;
        f)
            INPUT_FMT="$OPTARG"
            ;;
        O)
            OUT_DIR="$OPTARG"
            ;;
        n)
            BAM_NAME="$OPTARG"
            ;;
        l)
            LOGFILE="$OPTARG"
            ;;
        a)
            ALIGNMENT_TYPE="$OPTARG"
            ;;
        e)
            END_TO_END_PRESETS="$OPTARG"
            ;;
        c)
            LOCAL_PRESETS="$OPTARG"
            ;;
        N)
            NON_DETERMINISTIC="$OPTARG"
            ;;
        5)
            TRIM5="$OPTARG"
            ;;
        3)
            TRIM3="$OPTARG"
            ;;
        I)
            MINFRAGLEN="$OPTARG"
            ;;
        X)
            MAXFRAGLEN="$OPTARG"
            ;;
        t)
            THREADS="$OPTARG"
            ;;
        A)
            ADDITIONAL="$OPTARG"
            ;;
        h)
            HELP
            ;;
        :)
            echo ""$OPTARG" requires an argument"
            ;;
        \?) #unrecognized option - show help
            echo "Invalid option "$OPTARG""
            HELP
            ;;
    esac
done
 
#DEBUG
#echo "After getopts, the options are $*"
#echo "M1 is "$M1""
#echo "M2 is "$M2""
#echo -e "UNPAIRED is "$UNPAIRED"\n"

#It is common practice to call the shift command 
#at the end of your processing loop to remove 
#options that have already been handled from $@.
shift $((OPTIND -1))

#check for centrifuge image
if [[ ! -e "$SING_IMG" ]]; then
    echo "Missing SING_IMG \"$SING_IMG\""
    exit 1
fi

#building the fing argument string

#mutually exclusive
if [[ -n "$GENOME_DIR" ]]; then
    OPTSTRING="-g $GENOME_DIR"
elif [[ -n "$BT2_IDX" ]]; then
    OPTSTRING="-x $BT2_IDX"
else
    echo "Must have GENOME_DIR or BT2_IDX"
    exit 1
fi

if [[ -n "$ALIGNMENT_TYPE" ]] && [[ -n "$END_TO_END_PRESETS" ]]; then
    OPTSTRING="$OPTSTRING -a $ALIGNMENT_TYPE -e $END_TO_END_PRESETS"
elif [[ -n "$ALIGNMENT_TYPE" ]] && [[ -n "$LOCAL_PRESETS" ]]; then
    OPTSTRING="$OPTSTRING -a $ALIGNMENT_TYPE -L $LOCAL_PRESETS"
else
    echo "Must specify ALIGNMENT_TYPE and END_TO_END_PRESETS or LOCAL_PRESETS"
    exit 1
fi

#boolean
if [[ "$NON_DETERMINISTIC" -eq 1 ]]; then
    OPTSTRING="$OPTSTRING -N"
fi

#completely optional
build_opt_string -f $INPUT_FMT
build_opt_string -O $OUT_DIR
build_opt_string -t $THREADS
build_opt_string -l $LOGFILE
build_opt_string -5 $TRIM5
build_opt_string -3 $TRIM3
build_opt_string -I $MINFRAGLEN
build_opt_string -X $MAXFRAGLEN
build_opt_string -A $ADDITIONAL

#echo -e "After building the option string, we have:\n"
#echo -e ""$OPTSTRING"\n"

#Run bowtie

if [[ -z "$UNPAIRED" ]] && [[ -n "$M1" ]]; then

    echo -e "Looks like its just read pairs, no unpaired \n"
    IFS=' ' read -r -a M1ARRAY <<< "$M1"
    IFS=' ' read -r -a M2ARRAY <<< "$M2"

#    set -x
    for INDEX in "${!M1ARRAY[@]}"; do

        echo -e "Doing ${M1ARRAY[INDEX]}"
        echo -e "and ${M2ARRAY[INDEX]}\n"
        BAM_NAME=$(basename ${M1ARRAY[INDEX]} $INPUT_FMT).bam
        echo -e "Bam name will be $BAM_NAME\n"
        LOGFILE=$(basename ${M1ARRAY[INDEX]} $INPUT_FMT).log

        #this is where we would echo the command to a text file
        #that paramrun would then launch
        singularity exec $SING_IMG patric_bowtie2.py \
            -1 ${M1ARRAY[INDEX]} -2 ${M2ARRAY[INDEX]} \
            -l $LOGFILE -n $BAM_NAME $OPTSTRING

    done

elif [[ -n "$UNPAIRED" ]] && [[ -z "$M1" ]]; then

    IFS=' ' read -r -a UARRAY <<< "$UNPAIRED"

    for INDEX in "${!UARRAY[@]}"; do

        echo -e "Doing ${UARRAY[INDEX]}\n"
        BAM_NAME=$(basename ${UARRAY[INDEX]} $INPUT_FMT).bam
        echo -e "Bam name will be $BAM_NAME\n"
        LOGFILE=$(basename ${UARRAY[INDEX]} $INPUT_FMT).log

        singularity exec $SING_IMG patric_bowtie2.py \
          -U ${UARRAY[INDEX]} \
          -l $LOGFILE -n $BAM_NAME $OPTSTRING
          
    done

elif [[ -n "$UNPAIRED" ]] && [[ -n "$M1" ]]; then

    IFS=' ' read -r -a M1ARRAY <<< "$M1"
    IFS=' ' read -r -a M2ARRAY <<< "$M2"
    IFS=' ' read -r -a UARRAY <<< "$UNPAIRED"

    for INDEX in "${!UARRAY[@]}"; do

        echo -e "Doing ${M1ARRAY[INDEX]}"
        echo -e "and ${M2ARRAY[INDEX]}\n"
        echo -e "and ${UARRAY[INDEX]}\n"
        BAM_NAME=$(basename ${M1ARRAY[INDEX]} $INPUT_FMT).bam
        echo -e "Bam name will be $BAM_NAME\n"
        LOGFILE=$(basename ${M1ARRAY[INDEX]} $INPUT_FMT).log

        singularity exec $SING_IMG patric_bowtie2.py \
          -1 ${M1ARRAY[INDEX]} -2 ${M2ARRAY[INDEX]} \
          -U ${UARRAY[INDEX]} \
          -l $LOGFILE -n $BAM_NAME $OPTSTRING

    done

else

    echo "Something is wrong with how the reads were input"
    exit 1

fi

echo "Done with $0"
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

#!/bin/bash

#SBATCH -J bowtie2
#SBATCH -A iPlant-Collabs 
#SBATCH -N 4
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal

# Author: Scott G. Daniel <scottdaniel@email.arizona.edu>

###Uncomment when back on tacc#
#module load tacc-singularity 
#module load launcher
#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    
    echo "BAMS="./bowtie2-run.bam" #-b | --bams"
    echo "GFF_DIR="./genomes" #-g | --gff-dir"
    echo "GFF="./genomes/genome.gff" #-G | --gff"
    echo "LOG_FN="cuffdiff.log" #-l | --log-file"
    echo "THREADS="1" #-t | --threads"
    echo "OUT_DIR="./out_dir" #-O | --out-dir"
    echo "RRNA_GFF="./rRNA.gff" #-M | --mask-gff"
    echo "DEBUG="FALSE" #-d | --debug"
    echo "SILENT="FALSE" #-s | --silent"
    echo "MORE_ARGS="" #-A | --additional"

    echo "See cufflinks.py or make_graphs.R for additional help"
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

set -u

#
# Set up defaults for inputs, constants
#
BAMS="./bowtie2-run.bam" #-b | --bams"
GFF_DIR="./genomes" #-g | --gff-dir"
GFF="./genomes/genome.gff" #-G | --gff"
LOG_FN="cuffdiff.log" #-l | --log-file"
THREADS="1" #-t | --threads"
SING_IMG="cuffKeggR.img" 
OUT_DIR="./out_dir" #-O | --out-dir"
RRNA_GFF="./rRNA.gff" #-M | --mask-gff"
DEBUG="FALSE" #-d | --debug"
SILENT="FALSE" #-s | --silent"
MORE_ARGS="" #-A | --additional"

#Read the arguments
# In case you wanted to check what variables were passed
echo "ARG = $*"

while getopts :b:g:G:l:t:S:O:M:d:s:A:h ARG; do
    case $ARG in
        b)
            BAMS="$OPTARG"
            ;;
        g)
            GFF_DIR="$OPTARG"
            ;;
        G)
            GFF="$OPTARG"
            ;;
        l)
            LOG_FN="$OPTARG"
            ;;
        t)
            THREADS="$OPTARG"
            ;;
        o)
            OUT_DIR="$OPTARG"
            ;;
        M)
            RRNA_GFF="$OPTARG"
            ;;
        d)
            DEBUG=1
            ;;
        s)
            SILENT=1
            ;;
        A)
            MORE_ARGS="$OPTARG"
            ;;
        h)
            HELP
            ;;
        :)
            echo ""$OPTARG" requires an argument"
            ;;
        \?) #unrecognized option - show help
            echo "Invalid option "$OPTARG""
            HELP
            ;;
    esac
done

echo "What it looks like after parsing"
echo "ARG = $*"
#It is common practice to call the shift command 
#at the end of your processing loop to remove 
#options that have already been handled from $@.
shift $((OPTIND -1))


#If you have your own launcher setup on stampede2 just point MY_PARAMRUN at it
#this will override the TACC_LAUNCHER...
#PARAMRUN="${MY_PARAMRUN:-$TACC_LAUNCHER_DIR/paramrun}"

#check for centrifuge image
if [[ ! -e "$SING_IMG" ]]; then
    echo "Missing SING_IMG \"$SING_IMG\""
    exit 1
fi

#
# Verify existence of various directories, files
#
#[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"
# python script will do these checks

###Uncomment when back on TACC
#if [[ ! -d "$TACC_LAUNCHER_DIR" ]]; then
#    echo "Cannot find TACC_LAUNCHER_DIR \"$TACC_LAUNCHER_DIR\""
#    exit 1
#fi
#
#if [[ ! -f "$PARAMRUN" ]]; then
#    echo "Cannot find PARAMRUN \"$PARAM_RUN\""
#    exit 1
#fi


#Run cuffquant, cuffnorm, and cuffdiff
singularity run $SING_IMG \
-g $GFF_DIR \
-G $GFF \
-l $LOG_FN \
-t $THREADS \
-O $OUT_DIR \
-M $RRNA_GFF \
-d $DEBUG \
-s $SILENT \
-A $MORE_ARGS \
$BAMS

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

