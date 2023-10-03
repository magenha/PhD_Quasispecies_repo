import os

#Stablish the working directory
WD = os.system("WD=$pwd")
os.system("echo $WD")
#print(WD)
#Generate a report with fastqc for R1 reads
#bashCommand = f"find ${WD}/data -name '*R1*' "
#os.system(bashCommand)

#bashCommand = "mkdir test | touch test.txt"

#os.system(bashCommand)