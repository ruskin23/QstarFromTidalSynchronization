#!/bin/bash


python3 create_datafiles.py
cat good_data_file_* >GoodDataFile.txt
cat bad_data_file_* >BadDataFile.txt
rm good_data_file_* bad_data_file_*
