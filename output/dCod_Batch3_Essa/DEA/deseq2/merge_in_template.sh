#!/bin/bash

# get the complete gene list
## extract the first column of 3 fields to merge together
cut -f 1 -d ',' dea_station_4_station_1.csv > gene_id_list.txt
cut -f 1 -d ',' dea_station_4_station_2.csv >> gene_id_list.txt
cut -f 1 -d ',' dea_station_4_station_3.csv >> gene_id_list.txt

## sort and uniq
sort -u gene_id_list.txt > gene_id_sort.txt

# merge 3 files together with the sorted gene list based on the gene name
awk -F ',' 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' dea_station_4_station_1.csv gene_id_sort.txt | cut -f 2- -d ',' > merge_station_1.csv
awk -F ',' 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' dea_station_4_station_2.csv gene_id_sort.txt | cut -f 2- -d ',' > merge_station_2.csv
awk -F ',' 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' dea_station_4_station_3.csv gene_id_sort.txt | cut -f 2- -d ',' > merge_station_3.csv

# Put all 3 files above together with the list into one file
paste -d',' gene_id_sort.txt merge_station_1.csv merge_station_2.csv merge_station_3.csv > merge_stations.csv

