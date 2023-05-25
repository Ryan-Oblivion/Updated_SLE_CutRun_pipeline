# Updated_SLE_pipeline
use a directory of pe reads for this pipeline to analyze cut&amp;run data 

Once the full pipeline is completed, you will find the multiqc files in the rce directory.

This pipeline takes pe data that is in a directory but you will have to alter the glob pattern in the run_my_pe_pipeline.sh file.
You can do this by changing the parameter in the nextflow run command line. I will put a comment in the file to let you know where to do that.

Next, this pipeline does not download a human reference genome for you.
So once you download a reference genome go into the pipeline files and change all variables that contain a location to where my reference genome was, to where
yours is. I will see if i can update this later
