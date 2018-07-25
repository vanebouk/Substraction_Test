#!/bin/bash


#---------------------------------------------------
####################################################
# 				Assemblage
#---------------------------------------------------
####################################################

# files
#-------------------------------
pth=$(pwd)
pathfile="${pth}/file/"
echo $pathfile
pathID="./dataAssembly/"
echo ${pathID}

ls ./dataAssembly > ./file/Liste_SR.txt 

SR_List=$(cut -f1 -d'.' ${pathfile}Liste_SR.txt |sort | uniq) 
echo $SR_List

if [ -d "result_Assembly" ]
	then
		echo "result_Assembly existe !"
		rm -R result_Assembly
		mkdir result_Assembly
	else
		echo " result_Assembly no existe "
		mkdir result_Assembly
fi

if [ -d "data_Ass_SR" ]
	then
		echo "data_Ass_SR existe !"
		rm -R data_Ass_SR
		mkdir data_Ass_SR
	else
		echo " data_Ass_SR no existe "
		mkdir data_Ass_SR
	fi
#----------------
for Nom in $SR_List
	do
	echo ${Nom}
	echo -e "-----------------------------------------------------------\n">log.txt
	echo -e "Assemblage du fichier ${Nom} \n" >>log.txt
	echo -e "-----------------------------------------------------------\n">>log.txt
	echo -e "**** Assemblage ${Nom} **** \n">>log.txt
	dt=$('date')
	echo -e "Debut Assemblage ${Nom}: $dt \n" >>log.txt

	# recupérer les 3 premiers caractère du nom d'un fichier 
	id=$(echo ${Nom} | cut -c -2 )
	echo $id
	
	# kmers list
	kmersLength=("21" "33" "55")

	for kmers in "${kmersLength[@]}"
	do
		'date'
		dt=$('date')
		echo -e "	Debut Assemblage ${Nom} : $dt \n" >>log.txt
		#`rnaspades.py -k 21 --s1 ${pathID}${Nom}.fastq -o ${id}_Brut_k21_Contig_output`
		`rnaspades.py -k ${kmers} --s1 ${pathID}${Nom}.fastq -o ${Nom}_k${kmers}_Contig_output`
		
		dt=$('date')
		echo -e "	Fin Assemblage ${Nom} : $dt \n" >>log.txt
		
		# fasts to fastq
		#-----------------
		`mv ./${Nom}_k${kmers}_Contig_output/transcripts.fasta ./${Nom}_k${kmers}_Contig_output/${Nom}_k${kmers}.fasta`
		`perl ./fasta_to_fastq.pl ./${Nom}_k${kmers}_Contig_output/${Nom}_k${kmers}.fasta > ./data_Ass_SR/${Nom}_k${kmers}.fastq`
		# perl fasta_to_fastq.pl G6_SR.fasta > G6_SR.fastq
		

	done

	 `mv *_Contig_output result_Assembly`

done
