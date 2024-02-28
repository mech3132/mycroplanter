########### Read this for instructions on how to run process_athal_scan.py


First, make sure your image is CROPPED and that plant A1 is in the top right corner. Plant H12 should be in the bottom left.

Next, make sure all dependencies are installed. Run this to be sure:
pip install -r athalscan_requirements.txt

Finally, run the script. Be aware that python must invoke python3, not python2. I believe most computers these days come with default python3.

python process_athal_scan.py -i FILENAME.tif \
-m True \
-o output_directory

To see more options (e.g. to print the sub-images of each plant), run:
python process_athal_scan.py -h

 NOTE:
If progressbar doesn't install properly on mac, try running this in terminal:

conda install progressbar2
