#!/bin/bash

# # Start and end AA position are given as arguments. Translate to NT positions and calculate interval length
startAA=$1
endAA=$2
interval=$((3*($endAA - $startAA)))

newdirname="reconstruction_${startAA}_${endAA}"

# # Create directory to save reconstruction stuff if it does not exist
if [ -d "$newdirname" ]; then
echo "Network reconstruction has been performed."
exit
else
`mkdir -p $newdirname`;

echo "$newdirname directory has been created"

# # Copy files from Template folder to the new directory created
cp -r ./reconstruction_Template/. ./$newdirname/.

# # Enter new directory
cd $newdirname
# # Create Data directory
mkdir -p Data

# # # Move haplotypes file to Data from our haplotypes folder.
cd ../..
cp ./haplotypes/haplotypes_${startAA}_${endAA}/*.fasta ./reconstruction/${newdirname}/Data/.
# # # Enter Data folder
cd ./reconstruction/${newdirname}/Data

# # Copy names of the files wihout extension in fileNames.csv
for FILE in *; do echo ${FILE%.*} >> fileNames.csv; done
# # Exit to parent folder
cd ../

# # Run python scripts (Lui's code)
python3 script_findIDs.py
python3 script_extractAbundanceID.py

# Create metadata.csv file in Data folder
cd Data
cat nodesID.csv | wc -l >> metadata.csv
echo "0" >> metadata.csv
echo "$interval" >> metadata.csv
echo "" >> metadata.csv
echo "## metadata.csv" >> metadata.csv
echo "# nGenotypes in nodesID.csv." >> metadata.csv
echo "# refID." >> metadata.csv
echo "# seqLength." >> metadata.csv
cd ../

# # Run cpp script to generate network
./cScript_generateNetwork_

# Run python script to generate SubNetworks for each clade
python3 script_generateSubNetworks.py

fi
