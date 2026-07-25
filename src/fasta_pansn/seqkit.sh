#!/bin/bash

input=$1
output=$2

mkdir -p "$output"

for f in "$input"/*.fa.gz; do
    name=$(basename "$f" .fa.gz)
    seqkit replace -p '^>(.*)' -r "${name}#1#\${1}" "$f" > "${output}/${name}_pansn.fa"
done