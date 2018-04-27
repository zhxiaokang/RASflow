#!/bin/bash

# combine the countmerged.idx files based on the group that the samples belonging to

<<COMMENT
The group information is here:

control group (DMSO) (8 samples):
1, 6, 12, 18, 25, 31, 38, 44

BaP_10nM (7 samples):
2, 7, 13, 19, 26, 32, 45

BaP_1000nM (7 samples):
3, 8, 14, 20, 27, 33, 39

EE2_10nM (7 samples):
4, 9, 15, 21, 34, 40, 47

EE2_1000nM (6 samples):
5, 16, 22, 28, 35, 48

BaP_10nM + EE2_10nM (6 samples):
10, 23, 29, 36, 42, 49

BaP_1000nM + EE2_1000nM (6 samples):
11, 17, 24, 30, 37, 50
COMMENT

echo -e 'Gene_id\tsample_1\tsample_6\tsample_12\tsample_18\tsample_25\tsample_31\tsample_38\tsample_44' > control.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 1.countmerged.idx) <(awk -F'"' '{print $6}' 6.countmerged.idx) <(awk -F'"' '{print $6}' 12.countmerged.idx) <(awk -F'"' '{print $6}' 18.countmerged.idx) <(awk -F'"' '{print $6}' 25.countmerged.idx) <(awk -F'"' '{print $6}' 31.countmerged.idx) <(awk -F'"' '{print $6}' 38.countmerged.idx) <(awk -F'"' '{print $6}' 44.countmerged.idx) >> control.countmerged.idx

echo -e 'Gene_id\tsample_2\tsample_7\tsample_13\tsample_19\tsample_26\tsample_32\tsample_45' > BaP_low.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 2.countmerged.idx) <(awk -F'"' '{print $6}' 7.countmerged.idx) <(awk -F'"' '{print $6}' 13.countmerged.idx) <(awk -F'"' '{print $6}' 19.countmerged.idx) <(awk -F'"' '{print $6}' 26.countmerged.idx) <(awk -F'"' '{print $6}' 32.countmerged.idx) <(awk -F'"' '{print $6}' 45.countmerged.idx) >> BaP_low.countmerged.idx

echo -e 'Gene_id\tsample_3\tsample_8\tsample_14\tsample_20\tsample_27\tsample_33\tsample_39' > BaP_high.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 3.countmerged.idx) <(awk -F'"' '{print $6}' 8.countmerged.idx) <(awk -F'"' '{print $6}' 14.countmerged.idx) <(awk -F'"' '{print $6}' 20.countmerged.idx) <(awk -F'"' '{print $6}' 27.countmerged.idx) <(awk -F'"' '{print $6}' 33.countmerged.idx) <(awk -F'"' '{print $6}' 39.countmerged.idx) >> BaP_high.countmerged.idx

echo -e 'Gene_id\tsample_4\tsample_9\tsample_15\tsample_21\tsample_34\tsample_40\tsample_47' > EE2_low.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 4.countmerged.idx) <(awk -F'"' '{print $6}' 9.countmerged.idx) <(awk -F'"' '{print $6}' 15.countmerged.idx) <(awk -F'"' '{print $6}' 21.countmerged.idx) <(awk -F'"' '{print $6}' 34.countmerged.idx) <(awk -F'"' '{print $6}' 40.countmerged.idx) <(awk -F'"' '{print $6}' 47.countmerged.idx) >> EE2_low.countmerged.idx

echo -e 'Gene_id\tsample_5\tsample_16\tsample_22\tsample_28\tsample_35\tsample_48' > EE2_high.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 5.countmerged.idx) <(awk -F'"' '{print $6}' 16.countmerged.idx) <(awk -F'"' '{print $6}' 22.countmerged.idx) <(awk -F'"' '{print $6}' 28.countmerged.idx) <(awk -F'"' '{print $6}' 35.countmerged.idx) <(awk -F'"' '{print $6}' 48.countmerged.idx) >> EE2_high.countmerged.idx

echo -e 'Gene_id\tsample_10\tsample_23\tsample_29\tsample_36\tsample_42\tsample_49' > Mix_low.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 10.countmerged.idx) <(awk -F'"' '{print $6}' 23.countmerged.idx) <(awk -F'"' '{print $6}' 29.countmerged.idx) <(awk -F'"' '{print $6}' 36.countmerged.idx) <(awk -F'"' '{print $6}' 42.countmerged.idx) <(awk -F'"' '{print $6}' 49.countmerged.idx) >> Mix_low.countmerged.idx

echo -e 'Gene_id\tsample_11\tsample_17\tsample_24\tsample_30\tsample_37\tsample_50' > Mix_high.countmerged.idx
paste <(awk -F'"' '{print $2"\t"$6}' 11.countmerged.idx) <(awk -F'"' '{print $6}' 17.countmerged.idx) <(awk -F'"' '{print $6}' 24.countmerged.idx) <(awk -F'"' '{print $6}' 30.countmerged.idx) <(awk -F'"' '{print $6}' 37.countmerged.idx) <(awk -F'"' '{print $6}' 50.countmerged.idx) >> Mix_high.countmerged.idx
