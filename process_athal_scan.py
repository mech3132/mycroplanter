#!/bin/bash python3

####### IMPORTS #########

import argparse
from PIL import Image
import numpy as np  #for expanding array
#import re  for splitting via variable
#from random import choices
#from random import sample
import matplotlib.pyplot as plt
import os  # to make directories
import colorsys # for converting to hsv?
import math # for log, rounding
import sklearn # for pca
#from progressbar import ProgressBar ### This doesn't work on MACs, so removing 
from time import sleep

########## Set up ##########
parser = argparse.ArgumentParser(description = "For conducting coverage-based rarefaction. Finds minimum coverage amongst samples and rarefies all samples to that depth")
parser.add_argument('-i', '--input', type=str
                    , help='File path to image to be imported')
parser.add_argument('-p', '--parameters_filter', type=str
                    , help='Green and red ratio filters ( (p_i + 1) / (sum(p)+1) ) for plant mask. Input is comma separated green_float,red_float. [default: 0.38,0.355]', default='0.38,0.355')
parser.add_argument('-s', '--save_subimages', type=bool
                    , help='Save sub-images? [default: False]', default=False)
parser.add_argument('-b', '--save_blackandwhite', type=bool
                    , help='Save black and white image of plant filter? [default: False]', default=False)
parser.add_argument('-m', '--save_masked', type=bool
                    , help='Save masked image after plant filter? [default: False]', default=False)
parser.add_argument('-r', '--rownames', type=str
                   , help="From top to bottom, the name of each row of subimages. Number of names will determine how many segments the original image will be split into. Default is 12x8, with A1 in top RIGHT corner."
                   , default = ['1','2','3','4','5','6','7','8','9','10','11','12'])
parser.add_argument('-c', '--colnames', type=str
                   , help="From left to right, the name of each column of subimages. Number of names will determine how many segments the original image will be split into. Default is 12x8, with A1 in top RIGHT corner."
                   , default = ['H', 'G', 'F', 'E', 'D', 'C', 'B','A']
)
parser.add_argument('-o', '--output', type=str
                    , help='Directory to save final table and images', default="RGB_processing")				

args = parser.parse_args()
inputFP = args.input
save_subimages = args.save_subimages
save_bw = args.save_blackandwhite
save_masked = args.save_masked
params = args.parameters_filter
green_thresh, red_thresh = [float(x) for x in params.split(',')]
rowNames = args.rownames
colNames = args.colnames
outputFP = args.output


#### TESTING
#inputFP='age_concentration_n2c3001_cropped.tif'
## green_threshold = 0.36
## revblue_threshold = 0.70
#outputFP='./'
#save_subimages = False

#colNames = ['H', 'G', 'F', 'E', 'D', 'C', 'B','A']
#rowNames = ['1','2','3','4','5','6','7','8','9','10','11','12']

####### CLASSES ############

class image():
	def __init__(self, inputFP):
		# Load in image
		# Create classes A1 through H12
		# set thresholds
		self.im = Image.open(inputFP)
		# Make sure there's an alpha value
		self.im.putalpha(255)
		self.imarray = np.array(self.im)
	
	def make_custom_mask(self, maskname, maskarray):
		# For visualization later
		setattr(self, str(maskname), maskarray)
	
	def array_to_image(self, attr):
		temp = getattr(self,attr)
		return(Image.fromarray(temp))
	
	def mask_to_image(self, maskarrayname):
		## set image mask
		tempmask = getattr(self, maskarrayname)
		tempmaskreshape = tempmask.reshape(-2,tempmask.shape[-2])
		imarray_zeros = np.zeros(self.imarray.shape, dtype='uint8')
		imarray_zeros[:,:,3] = 255
		imarray_zeros[tempmaskreshape,:] = [255,255,255,255]
		return(Image.fromarray(imarray_zeros))
	
	def print_image(self, img_name, save_name):
		temp = getattr(self, img_name)
		temp.save(save_name+'.png')
	
	def divide_image(self, imarray, colNames, rowNames):
		self.colNames = colNames
		self.rowNames = rowNames
		ncol=len(colNames)
		nrow=len(rowNames)
		rpix = imarray.shape[0]
		cpix = imarray.shape[1]
		# n pix per group
		rwidth = rpix/nrow
		cwidth = cpix/ncol
		# Get ranges
		cranges =[[int(round(i*cwidth)),int(round((i+1)*cwidth-1))] for i in range(ncol)]
		rranges =[[int(round(i*rwidth)),int(round((i+1)*rwidth-1))] for i in range(nrow)]
		allSubimageNames = []
		print('Dividing images')
		sleep(0.2)
#		pbar = ProgressBar()
#		for nr,r in pbar(enumerate(rranges)): # doesn't work on macs
		maxrranges = len(rranges)
		for nr,r in enumerate(rranges):
			sleep(0.2)
			print('\r'+str(round(nr/maxrranges*100,0)) + '% done', end='')
			for nc,c in enumerate(cranges):
				subimage = imarray[r[0]:r[1], c[0]:c[1],:]
				subimageName = str('subimage_'+colNames[nc])+str(rowNames[nr])
				setattr(self, subimageName, subimage)
				allSubimageNames.append(subimageName)
		print('\r ========== 100% done ==========')
		self.subimageNames = allSubimageNames

	def process_all_subimages(self, function, params):
		# In this version, the FUNCTION outputs the number, rather than summing here. Allows multi-column outputs
		# If no attribute with list of masks, set here
		if hasattr(self, 'list_pixel_data'):
			getattr(self, 'list_pixel_data').append(function.__name__)
		else:
			newlist = [function.__name__]
			setattr(self, 'list_pixel_data', newlist)
			# Create new pixel data object if doesn't exist
		if hasattr(self, 'pixeldata'):
			None
		else:
			allPixData = {k:{} for k in self.subimageNames}
			setattr(self, 'pixeldata', allPixData)
		print('Processing '+function.__name__)
		sleep(0.2)
#		pbar = ProgressBar()
#		for i,s in pbar(enumerate(self.subimageNames)): # doesn't work on macs
		maxProg = len(self.subimageNames)
		for i,s in enumerate(self.subimageNames):
			sleep(0.2)
			print('\r'+str(round(i/maxProg*100,0)) + '% done', end='')
			temp = getattr(self, s)
			pixel_dat = function(temp, params)
			getattr(self, 'pixeldata')[s].update(pixel_dat)
		print('\r ========== 100% done ==========')
	   
	def print_subimages(self,output): ############ FIX TO BE ABLE TO PRINT B&W
		try:
			os.makedirs(outputFP+'/subimages')
		except FileExistsError:
			# directory already exists
			 pass
		print('Printing subimages')
		sleep(0.2)
#		pbar = ProgressBar()
#		for s in pbar(self.subimageNames): # doesn't work on macs
		maxProg = len(self.subimageNames)
		curr = 0
		for s in self.subimageNames:
			sleep(0.2)
			curr+=1
			print('\r'+str(round(curr/maxProg*100,0)) + '% done', end='')
			tempsubimage = getattr(self,s)
			Image.fromarray(tempsubimage).save(output+'/subimages/'+s+'.png')
		print('\r ========== 100% done ==========')
	
	def print_pixeldata(self, output):
		global outputFP
		colnames = list(self.pixeldata[list(self.pixeldata.keys())[0]].keys())
		toPrint = 'subimage\t'+'\t'.join(colnames)+'\n'
		for r in self.subimageNames:
			toPrint += r+'\t'
			for c in colnames:
				toPrint += str(self.pixeldata[r][c])+'\t'
			toPrint = toPrint.strip()
			toPrint +='\n'
		toPrint = toPrint.strip()
		f = open(outputFP+'/'+output+'.txt', 'w')
		f.write(toPrint)
		f.close()
		
######## FUNCTIONS ###########
def plant_mask(imagearray, green_thresh, red_thresh):
	## Make mask
	arrayr = np.array(imagearray[:,:,[0]], dtype=float)
	arrayg = np.array(imagearray[:,:,[1]], dtype=float)
	arrayb = np.array(imagearray[:,:,[2]], dtype=float)
	allPix = arrayr + arrayg + arrayb
    # Get ratio
	gratio = (arrayg+1)/(allPix+1)
	rratio = (arrayr+1)/(allPix+1)
	greenmask = gratio>green_thresh
	redmask = rratio>red_thresh
	blackmask = allPix>0
	nonwhitemask = (greenmask|redmask) & blackmask
	return(nonwhitemask)
def multiply_masks(masks):
    currentImage = np.array([])
    for m in masks:
        if len(currentImage.shape)==1:
            currentImage = m
        else:
            currentImage = currentImage*m
    return(currentImage)
def unite_masks(masks):
    currentImage = np.array([])
    for m in masks:
        if len(currentImage.shape)==1:
            currentImage = m
        else:
            currentImage = currentImage|m
    return(currentImage)
def summarise_rgb_pixels(imagearray, func):
    # THIS IGNORES ALL INSTANCES WHERE ALPHA IS ZERO
    # OUTPUT is a dictionary
    arrayr = np.array(imagearray[:,:,[0]], dtype=float).ravel()
    arrayg = np.array(imagearray[:,:,[1]], dtype=float).ravel()
    arrayb = np.array(imagearray[:,:,[2]], dtype=float).ravel()
    arrayA = np.array(imagearray[:,:,[3]], dtype=float).ravel()
    rflat = (arrayr[arrayA!=0]).tolist()
    gflat = (arrayg[arrayA!=0]).tolist()
    bflat = (arrayb[arrayA!=0]).tolist()
    summarise_pix = [func(x) for x in [rflat, gflat, bflat]]
    pix_dictionary = dict(zip(['pix_r_'+func.__name__,'pix_g_'+func.__name__,'pix_b_'+func.__name__], summarise_pix))
    return(pix_dictionary)
def count_all_pixels(imagearray, name_of_data):
    arrayA = np.array(imagearray[:,:,[3]], dtype=float).ravel()
    n = sum(arrayA!=0)
    n_dictionary = dict(zip([name_of_data],[n]))
    return(n_dictionary)
def calculate_ratios(imagearray,func):
    # THIS IGNORES ALL INSTANCES WHERE ALPHA IS ZERO
    arrayr = np.array(imagearray[:,:,[0]], dtype=float).ravel()
    arrayg = np.array(imagearray[:,:,[1]], dtype=float).ravel()
    arrayb = np.array(imagearray[:,:,[2]], dtype=float).ravel()
    arrayA = np.array(imagearray[:,:,[3]], dtype=float).ravel()
    allPix = arrayr + arrayg + arrayb
    rratio = (arrayr[arrayA!=0] / allPix[arrayA!=0]).tolist()
    gratio = (arrayg[arrayA!=0] / allPix[arrayA!=0]).tolist()
    bratio = (arrayb[arrayA!=0] / allPix[arrayA!=0]).tolist()
    # Use function to calculate mean, quantile, etc
    ratios = [func(x) for x in [rratio,gratio,bratio]]
    # Make into dictionary
    ratio_dictionary = dict(zip(['ratio_r_'+func.__name__,'ratio_g_'+func.__name__,'ratio_b_'+func.__name__], ratios))
    return(ratio_dictionary)
def calculate_differences(imagearray, func):
    arrayr = np.array(imagearray[:,:,[0]], dtype=float)
    arrayg = np.array(imagearray[:,:,[1]], dtype=float)
    arrayb = np.array(imagearray[:,:,[2]], dtype=float)
    arrayA = np.array(imagearray[:,:,[3]], dtype=float)
    # find differences,make into array, remove alpha=0, and then into list.
    diffrb = (arrayr-arrayb).ravel()[arrayA.ravel()!=0].tolist()
    diffgb = (arrayg-arrayb).ravel()[arrayA.ravel()!=0].tolist()
    diffgr = (arrayg-arrayr).ravel()[arrayA.ravel()!=0].tolist()
    # Use function to calculate differences
    diffs = [func(x) for x in [diffrb, diffgb, diffgr]]
    # Make into dictionary
    ratio_dictionary = dict(zip(['diff_rb_'+func.__name__,'diff_gb_'+func.__name__,'diff_gr_'+func.__name__], diffs))
    return(ratio_dictionary)

########### func for pixel cals
def average(x):
#	a=sum(x)/len(x)
	if len(x)>0:
		a=sum(x)/len(x)
	else:
		a=0
	return(a)
    
def median(x):
	if len(x)>0:
		a=np.quantile(x,0.5)
	else:
		a=0
	return(a)
    
def q5(x):
	if len(x)>0:
		a=np.quantile(x,0.05)
	else:
		a=0
	return(a)
    
def q95(x):
	if len(x)>0:
		a=np.quantile(x,0.95)
	else:
		a=0
	return(a)
	
#def average(x):
#	if len(x)>0:
#		a=sum(x)/len(x)
#	else:
#		a=np.nan
#	return(a)
#    
#def median(x):
#	if len(x)>0:
#		a=np.quantile(x,0.5)
#	else:
#		a=np.nan
#	return(a)
#    
#def q5(x):
#	if len(x)>0:
#		a=np.quantile(x,0.05)
#	else:
#		a=np.nan
#	return(a)
#    
#def q95(x):
#	if len(x)>0:
#		a=np.quantile(x,0.95)
#	else:
#		a=np.nan
#	return(a)

######## RUN #############

# Create and load image
try:
   os.makedirs(outputFP)
except FileExistsError:
   # directory already exists
   pass

## Load in image ##
imobj = image(inputFP)

######### MASK PLANT #############
imobj.make_custom_mask("plant_bw",plant_mask(imobj.imarray, green_thresh=green_thresh, red_thresh=red_thresh))
imobj.make_custom_mask("plant_masked", multiply_masks([imobj.plant_bw, imobj.imarray]))
imobj.divide_image(imobj.plant_masked, colNames, rowNames)

######## Calculate pixel data #########

# With plant mask, count the number of pixels
imobj.process_all_subimages(count_all_pixels, "all_plant_pixels") # store in pixel data file

# with plant mask, get median (can also use average) pixel ratios
#imobj.process_all_subimages(calculate_ratios, average) # store in pixel data file
imobj.process_all_subimages(calculate_ratios, q5) # store in pixel data file
imobj.process_all_subimages(calculate_ratios, median) # store in pixel data file
imobj.process_all_subimages(calculate_ratios, q95) # store in pixel data file

# With plant mask, find the average pixel
imobj.process_all_subimages(summarise_rgb_pixels, median) # store in pixel data file

# With plant mask, calculate pixel differences
imobj.process_all_subimages(calculate_differences, median) # store in pixel data file


#### Print images ####
if save_subimages:
    imobj.print_subimages(outputFP)
if save_bw:
    plant_bw = imobj.mask_to_image('plant_bw')
    plant_bw.save(outputFP+'/image_bw'+'.png')
if save_masked:
    plant_masked = imobj.array_to_image('plant_masked')
    plant_masked.save(outputFP+'/image_masked'+'.png')   

#### Print data ####
imobj.print_pixeldata('/pixel_data')

# print log


