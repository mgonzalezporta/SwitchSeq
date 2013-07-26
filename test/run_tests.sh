#!/bin/bash

## test1: get_switch - 4+4 samples
lorem -t get_switch \
	-s hsa -e 66 \
	-d /Users/mar/Desktop/app/data_10 \
	-i ./data_test1/expdata.txt \
	-o ./html_test1 \
	-c1 3-6 -c2 7-10

# -d /Users/aliquota/Desktop/data \
## test2: get_switch - 15+15 samples