#!/bin/bash

awk -F',' 'BEGIN{OFS="\t"} NR > 1 {printf "chr1\t%s\t%s\t%s__%s\n", $3-1, $4-1, $1, $2}' genes.csv