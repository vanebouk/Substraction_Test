#!/bin/bash


# paths
#--------

pth=$(pwd)/
#echo ${pth}
#cd ${pth}/
Pathfile="${pth}/file/"
pathBanque="${pth}Banques/"
pathresult="${pth}/result/"

# 
dt=$('date')
echo -e "	START DB Indexation : $dt \n" >${pathresult}log_indexDB.txt

# List of subtraction databases
#-------------------------------

Liste_Banques=$(cut -f1 -d'.' ${Pathfile}Liste_Banques.txt | uniq) 
#echo $Liste_Banques

# Script to free the RAM memory
chmod 755 Purge.sh

# Indexation of the subtraction databases
#-----------------------------------------

	for Banque in $Liste_Banques
	do
		echo -e "-----------------------------------------------------------\n"
		echo -e " Indexation of the subtraction databases : ${Banque}"
		echo -e "-----------------------------------------------------------\n"
			
		`bwa index ${pathBanque}${Banque}.fasta`
		./script/Purge.sh
	done
	
dt=$('date')
echo -e "	END DB Indexation : $dt \n" >${pathresult}log_indexDB.txt