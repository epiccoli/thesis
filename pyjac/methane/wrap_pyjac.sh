#!/bin/bash

echo Preparing pyjac python wrapper for gri30.cti file

python -m pyjac --lang c --input gri30.cti 2> err_lib.txt

python -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
