#!bin/bash python3

allOK = True

try: 
	import argparse
except:
	allOK = False
	print("Could not import argparse")
	
try: 
	from PIL import Image
except:
	allOK = False
	print("Could not import PIL")
	
try: 
	import numpy as np  #for expanding array
except:
	allOK = False
	print("Could not import numpy")
	
try: 
	import os  # to make directories
except:
	allOK = False
	print("Could not import os")
	
try: 
	import colorsys # for converting to hsv?
except:
	allOK = False
	print("Could not import colorsys")
	
try: 
	import math # for log, rounding
except:
	allOK = False
	print("Could not import math")
	
try: 
	import re # regular expression functions
except:
	allOK = False
	print("Could not import re")
	
try: 
	import sys # regular expression functions
except:
	allOK = False
	print("Could not import sys")
	
try: 
	from time import sleep 
except:
	allOK = False
	print("Could not import time")
	
try: 
	import pandas
except:
	allOK = False
	print("Could not import pandas")
	
if allOK:
	print("All requirements imported properly!") 
