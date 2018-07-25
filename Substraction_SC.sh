#!/bin/bash

# paths
#-------------------------------
dt=$('date')
echo -e " START Dual substraction : $dt \n" >>${pathresult}Substraction_log.txt
pth=$(pwd)

# script
#-------------------------------
pathScript="${pth}/script/"

# files
#-------------------------------
pathfile="${pth}/file/"

# data (Samplecriptomes files)
#-------------------------------
# a changer
pathdata="${pth}/data_Ass_SR/" 

# subtraction databases
#-------------------------------
pathBanque="${pth}/Banques/"

# Results
#-------------------------------
pathresult="${pth}/result_SC/"

# Total list of subtraction databases and IDs database
#------------------------------------------------------
Banque0='WolbachiaDB'
id0=$(echo $Banque0|cut -c -3)
echo ${id0}

Banque1="Wo2bachiaDB2"
#Banque1="AedesDB"
id1=$(echo $Banque1|cut -c -3)

Banque2='Wo3bachiaDB3'
#Banque2='AnophelesDB'
id2=$(echo $Banque2|cut -c -3)

#Banque3='CulexDB'
Banque3='Wo3bachiaDB3'
id3=$(echo $Banque3|cut -c -3)

#Banque4='DrosophilaDB'
Banque4='Wo2bachiaDB2'
id4=$(echo $Banque4|cut -c -3)

#Banque5='PedApisDB'
Banque5='Wo2bachiaDB2'
id5=$(echo $Banque5|cut -c -3)

#Banque6='Lutzomyia_PhlebotomusDB'
Banque6='WolbachiaDB'
id6=$(echo $Banque6|cut -c -3)

#Banque7='Tribolium_BombyxDB'
Banque7='Wo2bachiaDB2'
id7=$(echo $Banque7|cut -c -3)

Banque8='Wo3bachiaDB3'
#Banque8='IxodeDB'
id8=$(echo $Banque8|cut -c -3)

Banque9='Wo3bachiaDB3'
#Banque9='GlossinaDB'
id9=$(echo $Banque9|cut -c -3)

Banque10='Wo3bachiaDB3'
#Banque10='AutresDB'
id10=$(echo $Banque10|cut -c -3)

Banque11='Wo3bachiaDB3'
#Banque11='PrimatesDB'
id11=$(echo $Banque11|cut -c -3)

#-------------------------------
if [ -d "./result_SC" ]
then
	echo "result_SC existe !"
	rm -R result_SC
	mkdir result_SC
else
	echo " result_SC no existe "
	mkdir result_SC
fi

if [ -d "./data_SRSC" ]
then
 	echo "data_SRSC existe !"
	rm -R data_SRSC
 	mkdir data_SRSC
else
	echo " data_SRSC no existe "
	mkdir data_SRSC
fi


#---------------------------------------------------
####################################################

# 				SUBSTRACTION FUNCTION 

#---------------------------------------------------
####################################################

## Param, => BWA parameters 
## Banque, => database name
## Sample, => sample name
## IDdb, => database ID
## IDName => ID sample(transcriptome name)

Susbtraction_Function ()
{
	echo -e "substraction by mapping \n"
	echo -e "\nHere are the entered parameters : $1 $2 $3 $4 $5 $6\n"
	
	# $Param=$(echo $1)
	# echo $Param
	# $2=$Banque
	# $3=$Sample
	# $4=$IDName
	# $5=$IDdb
	# $6=$IDSamp
	# $7=$pth

	echo -e "-----------------------------------------------------------\n"
	echo -e " Mapping ${Sample}.fastq one ${Banque} \n"
	echo -e "-----------------------------------------------------------\n"

	`bwa mem ${Param} ${pathBanque}${Banque}.fasta ${pathdata}${Sample}.fastq > ${IDName}${IDdb}-aln-se.sam` 	 
	# -a  single-end 
	# -k = 19 minimum seed length  # longueur de k-mer  si diminu -> moins spécifique si augmente -> plus spécifique
	# -A = 1  Matching score : score des sequence matché
	# -c =10000 : IDNamebre d'occurrence 	
	# -B = 4  penalty for a mismatch # pénalise les mismatch 
	# -O = 6 gap open penalty : ouverture # pénaliser les gaps
	# -E =1 gap extension penalty; a gap of size k cost {-O} + {-E}*k [1] # penalise les grand gaps
	# -L = 5  penalty for clipping : pénalité de coupure  #

	./script/Purge.sh
	#---------------- Convert sam to bam ----------------------------------#
	#`samtools view -bS ${IDName}${IDdb}-aln-se.sam > ${IDName}${IDdb}.bam`	
	`samtools view -bS ${IDName}${IDdb}-aln-se.sam > ${IDName}${IDdb}.bam`		
	`samtools sort ${IDName}${IDdb}.bam ${IDName}${IDdb}_sorted`		
	`samtools index ${IDName}${IDdb}_sorted.bam`		
	`mv -f ${IDName}${IDdb}_sorted.bam.bai ${IDName}${IDdb}_sorted.bai`	

	#--------------- Mapped/Unmapped-----------------------------------------#
	#`samtools view -h ${IDName}${IDdb}_sorted.bam > ${IDName}${IDdb}_sorted.sam`
	./script/Purge.sh
	`samtools flagstat ${IDName}${IDdb}_sorted.bam > ${IDName}${IDdb}_sorted_stat.txt`
	./script/Purge.sh
	`samtools view -b -f 4 ${IDName}${IDdb}_sorted.bam > ${IDName}${IDdb}_Unmapped.bam`
	./script/Purge.sh
	`samtools view -b -F 4 ${IDName}${IDdb}_sorted.bam > ${IDName}${IDdb}_Mapped.bam`
	./script/Purge.sh
	#---------------Statistiques des Mapped/Unmapped ----------------------#
	`samtools flagstat ${IDName}${IDdb}_Unmapped.bam > ${IDName}${IDdb}_Unmapped_stat.txt`
	./script/Purge.sh
	`samtools flagstat ${IDName}${IDdb}_Mapped.bam > ${IDName}${IDdb}_Mapped_stat.txt`
	./script/Purge.sh
	`samtools mpileup -u -f ${pathBanque}${Banque}.fasta ${IDName}${IDdb}_Mapped.bam | bcftools call -c -v >${IDName}${IDdb}_Mapped.vcf`

	#-------------Convertion .bam to .fastq----------------------#
	`samtools bam2fq -nO -s ${IDName}${IDdb}.fastq ${IDName}${IDdb}_Unmapped.bam`
	./script/Purge.sh

	#------- Arrangement of the file results in the corresponding directories ------------#
	`mkdir ${IDSamp}_${IDdb}_Result_Subst`

	rm *.sam 	
	`mv *.bam ${IDSamp}_${IDdb}_Result_Subst`
	`mv *_stat.txt ${IDSamp}_${IDdb}_Result_Subst`
	`mv *.vcf ${IDSamp}_${IDdb}_Result_Subst`
	`mv *.bai ${IDSamp}_${IDdb}_Result_Subst`
	#rm *.bam 
	`cp ${IDName}${IDdb}.fastq ${IDSamp}_${IDdb}_Result_Subst`

}	


#---------------------------------------------------
####################################################

# 				APPLICATION

#---------------------------------------------------
####################################################
ls ./data_Ass_SR > ./file/Liste_SC.txt 

Liste_Sample=$(cut -f1 -d'.' ${pathfile}Liste_SC.txt |sort | uniq) 
echo $Liste_Sample
# $Sample => G6_SR_k21.fasta

# database number
Nbrdatabase=12

# boule sur le nombre d'échantillons
for Sample in $Liste_Sample
do	
	dt=$('date')
	echo -e " START substraction ${Sample}: $dt \n" >>${pathresult}Substraction_log.txt
	IDSamp=$(echo ${Sample} | cut -c -2 )
	
	pathdata="${pth}/data_Ass_SR/" 
	Sample0=$Sample

		#---#
		# 0 #
		#---#
		# parameters
		#------------
		Param="-a -k 14"
		Banque=${Banque0}
		IDdb=${id0}
		IDdb0=${id0}
		IDName=$(echo ${Sample} | cut -c -2 )
		IDName0=${IDName}
		echo $IDName0
		# call function
		#-----------------
		Susbtraction_Function "$Param" "$Banque" "$Sample" "id${i}" "$IDName" "$IDSamp"
		#---- out example ----#
		# Sample = G6_no_Ndup
		# IDName = G6
		# IDdb = Wol
		# G6WolUn.fastq

		#---#
		# 1 #
		#---#
		# parameters
		#------------
		Banque=${Banque1}
		IDdb=${id1}
		IDdb1=${IDdb}
		Sample="${IDName0}${IDdb0}"
		IDName="${Sample}"
		IDName1=${IDName}
		echo $IDName
		pathdata="./" 

		# call 
		#-------------------
		Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"

		#---#
		# 2 #
		#---#
		# parameters
		#------------
		Banque=${Banque2}
		IDdb=${id2}
		IDdb2=${IDdb}
		Sample="${IDName1}${IDdb1}" # G6WolWo2.fastq
		IDName="${Sample}"
		IDName2=${IDName}
		echo $IDName
		pathdata="./" 
		Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"

		# #---#
		# # 3 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque3}
		# IDdb=${id3}
		# IDdb3=${IDdb}
		# Sample="${IDName2}${IDdb2}" 
		# IDName="${Sample}"
		# IDName3=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"
				
		# #---#
		# # 4 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque4}
		# IDdb=${id4}
		# IDdb4=${IDdb}
		# Sample="${IDName3}${IDdb3}" 
		# IDName="${Sample}"
		# IDName4=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"


		# #---#
		# # 4 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque4}
		# IDdb=${id4}
		# IDdb4=${IDdb}
		# Sample="${IDName3}${IDdb3}" 
		# IDName="${Sample}"
		# IDName4=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"
	
		# #---#
		# # 5 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque5}
		# IDdb=${id5}
		# IDdb5=${IDdb}
		# Sample="${IDName4}${IDdb4}" 
		# IDName="${Sample}"
		# IDName5=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"
		
		# #---#
		# # 6 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque6}
		# IDdb=${id6}
		# IDdb6=${IDdb}
		# Sample="${IDName5}${IDdb5}" 
		# IDName="${Sample}"
		# IDName6=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"

		# #---#
		# # 7 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque7}
		# IDdb=${id7}
		# IDdb7=${IDdb}
		# Sample="${IDName6}${IDdb6}" 
		# IDName="${Sample}"
		# IDName7=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"

		# #---#
		# # 8 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque8}
		# IDdb=${id8}
		# IDdb8=${IDdb}
		# Sample="${IDName7}${IDdb7}" 
		# IDName="${Sample}"
		# IDName8=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"

		# #---#
		# # 9 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque9}
		# IDdb=${id9}
		# IDdb9=${IDdb}
		# Sample="${IDName8}${IDdb8}" 
		# IDName="${Sample}"
		# IDName9=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"
		
		# #---#
		# # 10 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque10}
		# IDdb=${id10}
		# IDdb10=${IDdb}
		# Sample="${IDName9}${IDdb9}" 
		# IDName="${Sample}"
		# IDName10=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"
		
		# #---#
		# # 11 #
		# #---#
		# # parameters
		# #------------
		# Banque=${Banque11}
		# IDdb=${id11}
		# IDdb11=${IDdb}
		# Sample="${IDName10}${IDdb10}" 
		# IDName="${Sample}"
		# IDName11=${IDName}
		# echo $IDName
		# pathdata="./" 
		# Susbtraction_Function "$Param" "$Banque" "$Sample" "$IDdb" "$IDName" "$IDSamp"


		if [ -d "${Sample0}_SRSC_Result" ]
		then
			echo "${Sample0}_SRSC_Result existe !"
			rm -R ${Sample0}_SRSC_Result
			mkdir ${Sample0}_SRSC_Result
		else
			echo " ${Sample0}_SRSC_Result no existe "
			mkdir ${Sample0}_SRSC_Result
		fi

		`mv *_Result_Subst ${Sample0}_SRSC_Result`
		`mv ${Sample0}_SRSC_Result result_SC`
		
		# rename files
		#---------------

		mv ${IDName2}${IDdb2}.fastq ${Sample0}_SRSC.fastq
		#`sed -e 's/^/${IDName2}${IDdb2}.fastq/${IDSamp}_SR.fastq/g'`
		#`rename -n 's/^/${IDName2}${IDdb2}.fastq/${IDSamp}_${IDdb2}_SR.fastq'`
		`cp ${Sample0}_SRSC.fastq data_SRSC`

		`rm *.fastq`

		dt=$('date')
		echo -e "	END Substraction ${Sample}: $dt \n" >>Substraction_log.txt
		'date'	
	done

#

dt=$('date')
echo -e "	End All Dual soustraction : $dt \n" >>${pathresult}Substraction_log.txt
'date'



#---------------------------------------------------
####################################################
# 				Contig Substraction
#---------------------------------------------------
####################################################










