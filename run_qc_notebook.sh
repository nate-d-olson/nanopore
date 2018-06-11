#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "You must specify output html file as argument to this script."
    exit
fi

jupyter nbconvert --execute --to html --TemplateExporter.exclude_input=True --TemplateExporter.exclude_output_prompt=True --TemplateExporter.exclude_input_prompt=True qc.ipynb --ExecutePreprocessor.timeout=600 --template full --output $1