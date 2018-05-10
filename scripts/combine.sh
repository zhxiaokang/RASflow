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
Reference site (station 4) (5 samples):
167512, 167912, 168112, 168412, 168612
Station 1 (4 samples):
160212, 160412, 160512, 160812
Station 2 (4 samples):
162712, 162812, 163112, 163512
Station 3 (3 samples):
164612, 165612, 166012
COMMENT

<<COMMENT
## A small Python code to help to generate the header for each output file
from __future__ import print_function  # if it's Python 2.7
list = [1, 2, 3]  # put the sample IDs for each group
group = 'control'  # fill in the group name
print('echo -e \'Gene_id', end = '')
for id in list:
    print('\\tsample_' + str(id), end = '')
print('\' > ' + group + '.countmerged.idx')
## To help to generate the command line of "paste"
print('paste <(awk -F\'\"\' \'{print $2}\' ' + str(list[0]) + '.countmerged.idx) ', end = '')
for id in list:
    print('<(awk -F\'\"\' \'{print $6}\' ' + str(id) + '.countmerged.idx) ', end = '')
print('>> ' + group + '.countmerged.idx\n')

COMMENT

echo -e 'Gene_id\tsample_1716412\tsample_1716512\tsample_1716812\tsample_1717012\tsample_1717812\tsample_1718012' > control.countmerged.idx
paste <(awk -F'"' '{print $2}' 1716412.countmerged.idx) <(awk -F'"' '{print $6}' 1716412.countmerged.idx) <(awk -F'"' '{print $6}' 1716512.countmerged.idx) <(awk -F'"' '{print $6}' 1716812.countmerged.idx) <(awk -F'"' '{print $6}' 1717012.countmerged.idx) <(awk -F'"' '{print $6}' 1717812.countmerged.idx) <(awk -F'"' '{print $6}' 1718012.countmerged.idx) >> control.countmerged.idx

echo -e 'Gene_id\tsample_1708012\tsample_1708612\tsample_1709112\tsample_1709212' > PAH_low.countmerged.idx
paste <(awk -F'"' '{print $2}' 1708012.countmerged.idx) <(awk -F'"' '{print $6}' 1708012.countmerged.idx) <(awk -F'"' '{print $6}' 1708612.countmerged.idx) <(awk -F'"' '{print $6}' 1709112.countmerged.idx) <(awk -F'"' '{print $6}' 1709212.countmerged.idx) >> PAH_low.countmerged.idx

echo -e 'Gene_id\tsample_1724712\tsample_1725412\tsample_1725712\tsample_1725912\tsample_1726812' > PAH_high.countmerged.idx
paste <(awk -F'"' '{print $2}' 1724712.countmerged.idx) <(awk -F'"' '{print $6}' 1724712.countmerged.idx) <(awk -F'"' '{print $6}' 1725412.countmerged.idx) <(awk -F'"' '{print $6}' 1725712.countmerged.idx) <(awk -F'"' '{print $6}' 1725912.countmerged.idx) <(awk -F'"' '{print $6}' 1726812.countmerged.idx) >> PAH_high.countmerged.idx

echo -e 'Gene_id\tsample_1709912\tsample_1710012\tsample_1710312\tsample_1710412\tsample_1710512' > PFAS_low.countmerged.idx
paste <(awk -F'"' '{print $2}' 1709912.countmerged.idx) <(awk -F'"' '{print $6}' 1709912.countmerged.idx) <(awk -F'"' '{print $6}' 1710012.countmerged.idx) <(awk -F'"' '{print $6}' 1710312.countmerged.idx) <(awk -F'"' '{print $6}' 1710412.countmerged.idx) <(awk -F'"' '{print $6}' 1710512.countmerged.idx) >> PFAS_low.countmerged.idx

echo -e 'Gene_id\tsample_1722912\tsample_1723512\tsample_1723812\tsample_1723912' > PFAS_high.countmerged.idx
paste <(awk -F'"' '{print $2}' 1722912.countmerged.idx) <(awk -F'"' '{print $6}' 1722912.countmerged.idx) <(awk -F'"' '{print $6}' 1723512.countmerged.idx) <(awk -F'"' '{print $6}' 1723812.countmerged.idx) <(awk -F'"' '{print $6}' 1723912.countmerged.idx) >> PFAS_high.countmerged.idx

echo -e 'Gene_id\tsample_167512\tsample_167912\tsample_168112\tsample_168412\tsample_168612' > station_4.countmerged.idx
paste <(awk -F'"' '{print $2}' 167512.countmerged.idx) <(awk -F'"' '{print $6}' 167512.countmerged.idx) <(awk -F'"' '{print $6}' 167912.countmerged.idx) <(awk -F'"' '{print $6}' 168112.countmerged.idx) <(awk -F'"' '{print $6}' 168412.countmerged.idx) <(awk -F'"' '{print $6}' 168612.countmerged.idx) >> station_4.countmerged.idx

echo -e 'Gene_id\tsample_160212\tsample_160412\tsample_160512\tsample_160812' > station_1.countmerged.idx
paste <(awk -F'"' '{print $2}' 160212.countmerged.idx) <(awk -F'"' '{print $6}' 160212.countmerged.idx) <(awk -F'"' '{print $6}' 160412.countmerged.idx) <(awk -F'"' '{print $6}' 160512.countmerged.idx) <(awk -F'"' '{print $6}' 160812.countmerged.idx) >> station_1.countmerged.idx

echo -e 'Gene_id\tsample_162712\tsample_162812\tsample_163112\tsample_163512' > station_2.countmerged.idx
paste <(awk -F'"' '{print $2}' 162712.countmerged.idx) <(awk -F'"' '{print $6}' 162712.countmerged.idx) <(awk -F'"' '{print $6}' 162812.countmerged.idx) <(awk -F'"' '{print $6}' 163112.countmerged.idx) <(awk -F'"' '{print $6}' 163512.countmerged.idx) >> station_2.countmerged.idx

echo -e 'Gene_id\tsample_164612\tsample_165612\tsample_166012' > station_3.countmerged.idx
paste <(awk -F'"' '{print $2}' 164612.countmerged.idx) <(awk -F'"' '{print $6}' 164612.countmerged.idx) <(awk -F'"' '{print $6}' 165612.countmerged.idx) <(awk -F'"' '{print $6}' 166012.countmerged.idx) >> station_3.countmerged.idx

