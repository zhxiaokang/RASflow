#!/bin/bash

# combine the countmerged.idx files based on the group that the samples belonging to

<<COMMENT
The group information is here:
control group for PAHs and PFAS exposure (6 samples):
1716412, 1716512, 1716812, 1717012, 1717812, 1718012
PAH_low (4 samples):
1708012, 1708612, 1709112, 1709212
PAH_high (5 samples):
1724712, 1725412, 1725712, 1725912, 1726812
PFAS_low (5 samples):
1709912, 1710012, 1710312, 1710412, 1710512
PFAS_high (4 samples):
1722912, 1723512, 1723812, 1723912
COMMENT

<<COMMENT
## A small Python code to help to generate the header for each output file
from __future__ import print_function  # if it's Python 2.7
list = [1716412, 1716512, 1716812, 1717012, 1717812, 1718012]  # fill in the sample IDs for each group
group = 'control'  # fill in the group name
print('echo -e \'Gene_id', end = '')
for id in list:
    print('\\tsample_' + str(id), end = '')
print('\' > ' + group + '.count')
## To help to generate the command line of "paste"
print('paste <(awk \'{print $1}\' ' + str(list[0]) + '.count) ', end = '')
for id in list:
    print('<(awk \'{print $3}\' ' + str(id) + '.count) ', end = '')
print('>> ' + group + '.count\n')

COMMENT

echo -e 'Gene_id\tsample_1716412\tsample_1716512\tsample_1716812\tsample_1717012\tsample_1717812\tsample_1718012' > control.count
paste <(awk '{print $1}' 1716412.count) <(awk '{print $3}' 1716412.count) <(awk '{print $3}' 1716512.count) <(awk '{print $3}' 1716812.count) <(awk '{print $3}' 1717012.count) <(awk '{print $3}' 1717812.count) <(awk '{print $3}' 1718012.count) >> control.count

echo -e 'Gene_id\tsample_1708012\tsample_1708612\tsample_1709112\tsample_1709212' > PAH_low.count
paste <(awk '{print $1}' 1708012.count) <(awk '{print $3}' 1708012.count) <(awk '{print $3}' 1708612.count) <(awk '{print $3}' 1709112.count) <(awk '{print $3}' 1709212.count) >> PAH_low.count

echo -e 'Gene_id\tsample_1724712\tsample_1725412\tsample_1725712\tsample_1725912\tsample_1726812' > PAH_high.count
paste <(awk '{print $1}' 1724712.count) <(awk '{print $3}' 1724712.count) <(awk '{print $3}' 1725412.count) <(awk '{print $3}' 1725712.count) <(awk '{print $3}' 1725912.count) <(awk '{print $3}' 1726812.count) >> PAH_high.count

echo -e 'Gene_id\tsample_1709912\tsample_1710012\tsample_1710312\tsample_1710412\tsample_1710512' > PFAS_low.count
paste <(awk '{print $1}' 1709912.count) <(awk '{print $3}' 1709912.count) <(awk '{print $3}' 1710012.count) <(awk '{print $3}' 1710312.count) <(awk '{print $3}' 1710412.count) <(awk '{print $3}' 1710512.count) >> PFAS_low.count

echo -e 'Gene_id\tsample_1722912\tsample_1723512\tsample_1723812\tsample_1723912' > PFAS_high.count
paste <(awk '{print $1}' 1722912.count) <(awk '{print $3}' 1722912.count) <(awk '{print $3}' 1723512.count) <(awk '{print $3}' 1723812.count) <(awk '{print $3}' 1723912.count) >> PFAS_high.count
