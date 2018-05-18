#!/usr/bin/env bash

#to be used to data munge those gff's

[[ $# -eq 0 ]] && echo "Usage: $0 [gff name]" && exit 1

DIR=$(dirname $1)
BASE=$(basename $1 .gff)
echo "Doing $DIR/$BASE.gff"

set -x
echo get rid of empty lines
egrep -v '^$' $1 > "$DIR/$BASE.temp"
echo get rid of comments
egrep -v '^#' $DIR/$BASE.temp > "$DIR/$BASE.temp2"
file "$DIR/$BASE.temp2"
echo only get CDS
grep -P "\tCDS\t" "$DIR/$BASE.temp2" > "$DIR/$BASE-CDS.gff"
echo "get the rRNA for cuffdiff (or other progs) that support filtering by rRNA seq"
grep -P "\trRNA\t" "$DIR/$BASE.temp2" > "$DIR/$BASE-rRNA.gff" 

echo removing temp files if everything went sucessfull
[[ -s $DIR/$BASE-CDS.gff ]] && rm $DIR/$BASE.temp
[[ -s $DIR/$BASE-CDS.gff ]] && rm $DIR/$BASE.temp2
