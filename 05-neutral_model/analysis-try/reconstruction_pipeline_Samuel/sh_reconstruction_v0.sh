#!/bin/bash

# # Start and end AA position are given as arguments. Translate to NT positions and calculate interval length
#protein=$1
#startAA=$2
#endAA=$3
#interval=$((3*(1+($endAA - $startAA))))

WD=`pwd`

newdirname="reconstruction_qbeta_simulated"

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
cd $WD
#echo "($pwd)"
cp ./Data/*.fasta ./${newdirname}/Data/.
# # # Enter Data folder
cd ./${newdirname}/Data

# # Copy names of the files wihout extension in fileNames.csv
for FILE in *; do echo ${FILE%.*} >> fileNames.csv; done
# # Exit to parent folder
cd ../

# exit 0 # # to exit the program
# kill -9 $$ # # to kill the process

echo "Executing script_findIDs.py"
# # Create nodesID.csv file, one line for each haplo seq identified.
python3 script_findIDs.py
# # Run python scripts (Iker's code) -
# # same result but we create a dict with all haplo identified by order
# # To be considered for temporal network reconstruction
# python3 script_findIDs_iker.py


echo "Executing script_extractAbundanceID.py"
# # Create abundances files for each clades
python3 script_extractAbundanceID.py

echo "Executing py_extractAbundanceGN.py"
# # Create GN abundance json file
python3 py_extractAbundanceGN.py #$protein $startAA $endAA

echo "Creating metadata.csv"
# Create metadata.csv file in Data folder
cd Data
cat nodesID.csv | wc -l >> metadata.csv
echo "0" >> metadata.csv
echo "201" >> metadata.csv
echo "" >> metadata.csv
echo "## metadata.csv" >> metadata.csv
echo "# nGenotypes in nodesID.csv." >> metadata.csv
echo "# refID." >> metadata.csv
echo "# seqLength." >> metadata.csv
cd ../


echo "Running cScript_generateNetwork_"

# # Run cpp script to generate network
g++-11 -std=c++17 -o cScript_generateNetwork_ cScript_generateNetwork_.cpp
./cScript_generateNetwork_

echo "Running script_generateSubNetworks.py"

# Run python script to generate SubNetworks for each clade
python3 script_generateSubNetworks.py

fi
