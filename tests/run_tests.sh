#!/bin/bash

## test1: get_switch - 4+4 samples
echo "Running test 1..."
switchseq -t get_switch \
	-s hsa -e 66 \
	-d ./data_annot \
	-i ./data_test1/expdata.txt \
	-o ./html_test1 \
	-c1 3-6 -c2 7-10

## test2: get_switch - 12+12 samples
echo "Running test 2..."
switchseq -t get_switch \
	-s hsa -e 66 \
	-d ./data_annot \
	-i ./data_test2/expdata.txt \
	-f ./data_test2/significant_dtu.txt \
	-o ./html_test2 \
	-c1 3-14 -c2 15-26

## test3: get_switch - 4+4 samples (with technical replicates)
echo "Running test 3..."
switchseq -t get_switch \
        -s hsa -e 66 \
        -d ./data_annot \
        -i ./data_test3/expdata.txt \
        -o ./html_test3 \
        -c1 3-7 -c2 8-12

## test4: get_switch - 2+2 samples; UCSC annotation
echo "Running test 4..."
switchseq -t get_switch \
        -s hsa -e 66 \
	-a TRUE \
        -d ./data_annot \
        -i ./data_test4/expdata.txt \
        -o ./html_test4 \
        -c1 3-4 -c2 5-6
