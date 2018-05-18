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

