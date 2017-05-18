#!/bin/bash

echo Preparing pyjac python wrapper for USC.cti file

python -m pyjac --lang c --input USC.cti 2> err_lib.txt

python -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
