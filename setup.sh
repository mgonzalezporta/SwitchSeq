#!/bin/bash

checkPerlModule() {
    perl -M$1 -e '' > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo "   $1 - OK"
    else
        echo "   $1 - NOT FOUND"
	echo "   Installing $1..."
        perl -MCPAN -e 'install $1'
    fi
}

checkBin() {
    if hash $1 2>/dev/null; then
        echo "   $1 - OK"
    else
        echo "   $1 - NOT FOUND"
	echo "   Please install $1 manually (see https://github.com/mgonzalezporta/lorem/wiki/Setup-guide)"
    fi
}

checkRPackage() {
result=`R --no-save --quiet <<RSCRIPT
is.element('$1', installed.packages()[,1])
RSCRIPT`

    if [[ $result == *TRUE* ]]; then
        echo "   $1 - OK"
    else
        echo "   $1 - NOT FOUND"
	echo "   Installing $1..."
	git clone https://github.com/mgonzalezporta/ipsum.git
	cd ipsum/
	R CMD build Ipsum
	R CMD install Ipsum_1.0.tar.gz
	rm Ipsum_1.0.tar.gz
    fi
}

echo "# Checking Perl modules..."
modules=( 'Getopt::Long' 'Pod::Usage' 'Exporter' 
	'LWP::Simple' 'List::Util' 'Text::Template' 'JSON')

for m in "${modules[@]}"
do
	checkPerlModule $m
done

echo "# Checking needle installation..."
checkBin needle

echo "# Checking R package Ipsum..."
checkRPackage Ipsum
