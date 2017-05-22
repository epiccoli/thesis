#!/bin/bash

echo Preparing pyjac python wrapper for h2o2mech_Lietal_2003.cti file

python -m pyjac --lang c --input h2o2mech_Lietal_2003.cti 2> err_lib.txt

python -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
