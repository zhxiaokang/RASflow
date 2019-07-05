#!/bin/bash

# format the count file as standard input for combine2group.py
# two columns: ID (\t) count
# file named as: sample_count.tsv

# two parameters: the original count file sample.countmerged.idx, and the output file

awk -F'"' '{print $2"\t"$6}' $1 > $2
