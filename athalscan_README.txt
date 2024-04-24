########### Read this for instructions on how to run process_athal_scan.py


First, make sure your image is CROPPED and that plant A1 is in the top right corner. Plant H12 should be in the bottom left.

Next, make sure all dependencies are installed. Run this to be sure:

python3 -m pip install -r athalscan_requirements.txt

Now, let's check your installations. Run:

python3 testing_installations.py

If you get a message saying all installations have been imported properly, you're good to go! If not, then try to manually install the missing packages by doing the following:

python3 -m pip install PACKAGENAME

You may need to google what kinds of solutions people have had if you are still running into issues.


Finally, run the script. To run the most basic version of the script, use the following code:

python3 process_athal_scan.py -i sample_image.tif -o output_directory



To see more options (e.g. to print the sub-images of each plant, set custom functions, etc), run:
python3 process_athal_scan.py -h

######

After running process_athal_scan.py, you may want to merge your metadata (which treatments are in which wells?) and your data (pixel data).

To merge, run the following (assuming your metadata is in tab-delimited format):


python3 merge_athal_and_meta.py -p output_directory/pixel_data.txt -m sample_metadata.tsv -o merged_data.txt

Again, for more options and to change the input/output formats, you can run:
python3 merge_athal_and_meta.py -h

