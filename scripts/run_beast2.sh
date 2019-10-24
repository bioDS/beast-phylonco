#!/bin/bash
# Description: runs beast MCMC with the phylonco package
# Usage: ./run_beast2.sh [options] <xml file>
file_path=${@:$#}
options=${*%${!#}}
current_path="$(pwd)"
jar_path="$current_path/../dist/phylonco.jar"
output_path="output/beast"

if ! [[ -d $output_path ]]; then mkdir $output_path; fi
if ! [[ $file_path = /* ]]; then file_path="$current_path/$file_path"; fi

cd $output_path
echo "writing output to: $output_path"

if [[ -d $file_path ]]; then
	for f in "$file_path/"*.xml; do
		echo "java -jar $jar_path $options $f"
		java -jar "$jar_path" "$options" "$f" 
	done	
else
	echo "java -jar $jar_path $options $file_path"
	java -jar "$jar_path" "$options" "$file_path"
fi

echo "done!"