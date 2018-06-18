This nanopore pipeline was developed for processing GIAB nanopore data. The following documentation is intended for group members running the pipeline on the Stanford Sherlock cluster. 

# Set-up
Create a python 3.6 virtualenv for installing pipeline dependencies. The following commands will work on sherlock:
```
## Load required modules - ONLY ON SHERLOCK
module load readline R python/3.6.1
## Create and activate virtualenv
virtualenv -p python3 nanopore_env
source nanopore_env/bin/activate
```

The python package poetry is used to install this package and its dependencies.

Run the following commands in the python 3.6 virtualenv:

```
curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
git clone https://github.com/nspies/nanopore.git
cd nanopore
poetry install
```



## Installing Depencies 

### Albacore basecaller.  
For server side installations, download and install the albacore basecaller from [here](https://community.nanoporetech.com/downloads). This is not necessary on the sequencing machines.
Link to Albacore wheel file is available on the Nanopore community page under software downloads. 
You will need to register for an account to access the community page. 

Update the following commands with links to the latest version of albacore:

```
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-<version>-cp36-cp36m-manylinux1_x86_64.whl
pip install ont_albacore-<version>-cp36-cp36m-manylinux1_x86_64.whl
```

### Minimap2
For mapping long ONT reads. 
```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```


## Accessing data
Contact the group data administrator to get access to OAK. 

With the dependencies installed and appropriate permissions you can now run the pipeline. 

# Sequencer-side

## Copying to server

To run the copying script: `python copy-to-server.py --user <user_name> /Library/MinKNOW/data/reads <sample_name> ~/temp`, where **sample_name** must be the name of the run specified in MinKNOW. For example, if your output files are ending up in the sub-directory `20180124_2338_180124_MyRunName`, you would specify `MyRunName`. All runs using that name will be uploaded, allowing for restarts (the same name must be used each subsequent re-run).

The script can be started before, during, or after a run. It will continue running until cancelled using ctrl-c. It will copy files over once enough have been produced (in batches of at least 1000), and will copy over remaining files after sequencing has stopped running. The `~/temp` directory is where the raw fast5 files will be moved to, then combined together into .tar files which are subsequently rsync'ed to the server.

# Server-side

## Running

All of the following require the activation of the virtualenv created above. 
The following modules are required for running the pipeline:

```
module load python/3.6.1 R biology samtools/1.8
```

__Note__ make sure to load python/3.6.1 module before activating the virtualenv. 

### Per-sample basecalling and alignment scripts

To run: `python main.py [flowcell-id]`. If no flowcell ID is specified, all samples will be run, skipping those that appear to have run to completion.

This requires the flowcell to have been entered properly into the google spreadsheet, and the data to have been transferred over to oak storage.

### Per-sample QC

To run: `bash run_qc_notebook.sh <qc.html>`, where you should replace **qc.html** with the appropriate output file for the formatted results notebook. Execute script within repo directory or have `qc.ipynb` in your path.  

### Cutting a release

To run: `bash cut_release.sh <release_name>`, where **release_name** should be something like "ultralong_combined_2018-06-11". This will combine fast5 files, fastq files and remap against hg19 and hg38 producing cram files.
