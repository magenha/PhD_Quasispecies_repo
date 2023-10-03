#!/bin/bash

# # Start and end AA position are given as arguments. Translate to NT positions and calculate interval length
#protein=$1
#startAA=$2
#endAA=$3
#interval=$((3*(1+($endAA - $startAA))))
WD=`pwd`

newdirname="reconstruction_qbeta_simulated"


cd $newdirname


echo "Creating metadata.csv"
# Create metadata.csv file in Data folder
cd Data
cat nodesID.csv | wc -l >> metadata.csv
echo "0" >> metadata.csv
echo "271" >> metadata.csv
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
