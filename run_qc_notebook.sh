#!/usr/bin/env bash

jupyter nbconvert --execute --to html --TemplateExporter.exclude_input=True --TemplateExporter.exclude_output_prompt=True --TemplateExporter.exclude_input_prompt=True qc.ipynb --template full --output $1