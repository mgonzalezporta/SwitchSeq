#!/bin/bash

#needle

checkPerlModule() {
    perl -M$1 -e '' > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "   $1 - OK"
    else
        echo "   $1 - NOT FOUND"
    fi
}

echo "# Checking Perl modules..."
modules=( 'Getopt::Long' 'Pod::Usage' 'Exporter' 
	'LWP::Simple' 'List::Util' 'Text::Template' 'JSON')

for m in "${modules[@]}"
do
	checkPerlModule $m
done