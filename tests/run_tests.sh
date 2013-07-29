#!/bin/bash

## test1: get_switch - 4+4 samples
echo "Running test 1..."
lorem -t get_switch \
	-s hsa -e 66 \
	-d /Users/aliquota/Desktop/data \
	-i ./data_test1/expdata.txt \
	-o ./html_test1 \
	-c1 3-6 -c2 7-10

# 	-d /Users/mar/Desktop/app/data_10 \
# -d /Users/aliquota/Desktop/data \

## test2: get_switch - 12+12 samples
echo "Running test 2..."
lorem -t get_switch \
	-s hsa -e 66 \
	-d /Users/aliquota/Desktop/data \
	-i ./data_test2/expdata.txt \
	-o ./html_test2 \
	-c1 3-14 -c2 15-26