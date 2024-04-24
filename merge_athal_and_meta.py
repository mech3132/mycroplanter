#!bin/bash

import pandas as pd
import numpy as np
import argparse
import os  # to make directories
import sys # for errors
from time import sleep

### Get inputs ####
parser = argparse.ArgumentParser(description = "Combines output from process_athal_scan.py and a csv (or tsv) metadata file. This combines tables by row and col ONLY, so you MUST have one column that is 'col' (Capitilized letters) and one column that is 'row' (numbers; no leading zeros) in your metadata. You may also merge multiple plates of data at once (multiple pixeldata files and multiple metadata files) BUT then you must also include 'plate' was one of your column headers. You can add a 'plate' column header using process_athal_scan.py with the add_variable argument.")
parser.add_argument('-p', '--pixeldata', type=str
                    , help='File path to process_athal_scan.py output. If multiple files, separate with commas.')
parser.add_argument('-m', '--meta', type=str
					, help="File path to metadata. Default format is tab-delimited. To change delimiter, use -d flag. If multiple files, separate with commas")
parser.add_argument('-d', '--delimiter', type=str
					, help = "Default INPUT delimiter is tab. For csv, input ','.", default='\t')
parser.add_argument('-D', '--output_delimiter', type=str
					, help = "Default OUTPUT delimiter is tab. For csv, input ','.", default='\t')
parser.add_argument('-o', '--output', type=str
                    , help='File path to save merged files')

args = parser.parse_args()
pixeldataFP = args.pixeldata.split(',')
metaFP = args.meta.split(',')
delim = args.delimiter
odelim = args.output_delimiter
output = args.output

## For testing
#metaFP = '00_rawdata/meta/metadata_plateA.tsv,00_rawdata/meta/metadata_plateB.tsv'.split(',')
#pixeldataFP = '01_processandmerge/plateA/pixel_data.txt,01_processandmerge/plateB/pixel_data.txt'.split(',')
#delim = ','
#odelim = '\t'
#output='merged_dat.txt'
############## Start #################

# iterate through all pixeldata; concatenate
dat = pd.DataFrame()
for pixdat in pixeldataFP:
	if dat.shape[0]<1:
		dat = pd.read_csv(pixdat, sep='\t')
	else:
		newdat = pd.read_csv(pixdat, sep='\t')
		sharedCol = [x for x in set([x for x in newdat.columns] + [x for x in + dat.columns])]
		try:
			dat = pd.merge(dat, newdat, on=sharedCol, how='outer')
		except: 
			print("Warning: there was not a plate column in pixel data, which means you cannot differentiate between plates that were merged together. Check that your outputs have a plate column labelled 'plate'")
		
	
meta = pd.DataFrame()
for met in metaFP:
	if meta.shape[0]<1:
		meta = pd.read_csv(met, sep='\t')
	else:
		newmeta = pd.read_csv(met, sep='\t')
		sharedColMet = [x for x in set([x for x in newmeta.columns] + [x for x in + meta.columns])]
		try:
			matchDtypes = meta.dtypes[sharedColMet] != newmeta.dtypes[sharedColMet]
			while matchDtypes.any(0):
				change=matchDtypes[matchDtypes==True].index[0]
				if newmeta.dtypes[change] == np.object_ or meta.dtypes[change] == np.object_ :
					newmeta = newmeta.astype({change: np.object_})
					meta = meta.astype({change: np.object_})
					print("WARNING: changed column "+change+" to a string")
					matchDtypes = meta.dtypes[sharedColMet] != newmeta.dtypes[sharedColMet]
				if newmeta.dtypes[change] == np.float64 or meta.dtypes[change] == np.float64 :
					newmeta = newmeta.astype({change: np.float64})
					meta = meta.astype({change: np.float64})
					print("WARNING: changed column "+change+" to a float")
					matchDtypes = meta.dtypes[sharedColMet] != newmeta.dtypes[sharedColMet]
				if newmeta.dtypes[change] != meta.dtypes[change]: # Catchall; change into string if unsure
					newmeta = newmeta.astype({change: np.object_})
					meta = meta.astype({change: np.object_})
					print("WARNING: changed column "+change+" to a string")
					matchDtypes = meta.dtypes[sharedColMet] != newmeta.dtypes[sharedColMet]
			meta = pd.merge(meta, newmeta, on=sharedColMet, how='outer')
		except: 
			print("Warning: there was not a plate column in metadata, which means you may not be able to differentiate between plates that were merged together. ")

# Check that dataframes are loaded properly
if meta.shape[1]<=1:
	print("\WARNING: metadata has only one column. Are you sure you used the correct delimiter?\n")
if dat.shape[1]<=1:
	print("\WARNING: pixel data has only one column. Have you confirmed it is a tab-delimited file?\n")

# Split subimage names
try:
	dat[['remove','well']] = dat['subimage'].str.split('_',expand=True)
	dat['row'] = dat['well'].str.slice(stop=1)
	dat['col'] = dat['well'].str.slice(start=1)
	dat = dat.drop(['remove','subimage','well'], axis=1)
except:
	print('\nCould not find subimage column in pixel data. Check that it still exists?\n')
	sys.exit(1)

# Change row and col types
dat[['row','col']] = dat[['row','col']].astype('str')
# Do the same to meta
try:
	meta[['row','col']] = meta[['row','col']].astype('str')
except:
	print("\nCould not convert the variables 'row' or 'col' to strings to merge tables. Check if they exist in your metadata?\n")
	sys.exit(1)

####### MERGE #########
try:
	alldat = pd.merge(meta, dat, on=['row','col','plate'], how='outer')
except:
	try:
		alldat = pd.merge(meta, dat, on=['row','col'], how='outer')
		if len(pixeldataFP)>1:
			print("Warning: more than one pixel data spreadsheet was provided but there was no 'plate' column to merge by. Please check output to see if merging went wrong")
	except:
		print("ERROR: row and col do not exist for merging. Please check that your metadata file includes these two variables, and that your pixel_data file includes the column 'subimage'.")
		sys.exit(1)

# Check merge
if alldat.shape[0] != meta.shape[0]:
	print("Warning: the number of metadata entries does not equal the number of final data points. Please ensure your metadata file contains all the samples you need, and that row and col values are correct. If you purposely provided a metadata table that is larger or smaller than your pixel data table, be aware that merging will KEEP all entries and that you should filter out NAs downstream.")


print("\n\nMerged data and metadata.")

alldat.to_csv(output, sep=odelim, index=False)
