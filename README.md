## Set-up

Run the following commands in a python 3.6 virtualenv:

```
curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
poetry install
```

Then install the albacore basecaller from [here](https://community.nanoporetech.com/downloads).

## Running

All of the following require the activation of the virtualenv created above.

### Per-sample basecalling and alignment scripts

To run: `python main.py [flowcell-id]`. If no flowcell ID is specified, all samples will be run, skipping those that appear to have run to completion.

This requires the flowcell to have been entered properly into the google spreadsheet, and the data to have been transferred over to oak storage.

### Per-sample QC

To run: `bash run_qc_notebook.sh <qc.html>`, where you should replace **qc.html** with the appropriate output file for the formatted results notebook. 

### Cutting a release

To run: `bash cut_release.sh <release_name>`, where **release_name** should be something like "ultralong_combined_2018-06-11". This will combine fast5 files, fastq files and remap against hg19 and hg38 producing cram files.