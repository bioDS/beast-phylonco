#!/bin/bash
# counts number of sites per read depth
# usage: ./count_coverage.sh <path>
# output: number of sites, depth
for file in $1
do
	echo "processing file: $file"
	cat $file | awk '{print $3}' | sort -n | uniq -c > "$file.count"
done