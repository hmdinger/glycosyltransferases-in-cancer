import os
import sys

filt = 'de_gts_master_filtered.tsv'
fc_filt = 'de_gts_master_fc_filtered.tsv'
fc_high = 'de_gts_master_fc_filtered_high.tsv'
fc_higher = 'de_gts_master_fc_filtered_higher.tsv'
fc_highest = 'de_gts_master_fc_filtered_highest.tsv'

filtered = [filt, fc_filt, fc_high, fc_higher, fc_highest]

for file in filtered:
	print(file)
	file = file.replace("de_gts_master_", "")
	print(file)
	file = file.replace(".tsv", "")
	print(file)
