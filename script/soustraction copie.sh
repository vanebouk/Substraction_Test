#!/bin/bash

dt=$('date')
echo -e "	Debut All soustraction k14 : $dt \n" >log.txt
pth=$(pwd)/
echo ${pth}
cd ${pth}/

Pathfile="${pth}/file/"
File="/Users/boukoubavanesa/Documents/Data_NB/Assemblage/files_assembly/"
file2="/Users/boukoubavanesa/Documents/Data_NB/Soustraction/file/"
pathBanque="/Users/boukoubavanesa/Documents/Data_NB/Soustraction/Banques/"
#/Users/boukoubavanesa/Documents/Data_NB/Soustraction/Banques/
pathScript="${pth}script/"

file3="/Users/boukoubavanesa/Documents/Data_NB/Soustraction/data/"
pathAssemblage="/Users/boukoubavanesa/Documents/Data_NB/Assemblage/"

# Test si les fichiers existe 

if [ ! -e ${Pathfile}Liste_Transcriptomes.txt ] 
then
	echo -e "Le fichier n'esxiste pas \nVeiellez verifier que le fichier est dans le meme repertoire que le script \nMERCI"
else

Liste_noms=$(cut -f1 -d'.' ${Pathfile}Liste_Transcriptomes.txt |sort | uniq) 
#Liste_Banques=$(cut -f1 -d'.' ${Pathfile}Liste_Banques.txt | uniq) 
Banque=$(cut -f1 -d'.' ${Pathfile}Banque.txt | uniq)
#Banque2=$(cut -f1 -d'.' Banque2.txt | uniq)

#=================Indexation des banques de soustraction================================
	
Banque0='WolbachiaDB'
Banque1="AedesDB"
Banque2='AnophelesDB'
Banque3='CulexDB'
Banque4='DrosophilaDB'
Banque5='PedApisDB'
Banque6='Lutzomyia_PhlebotomusDB'
Banque7='Tribolium_BombyxDB'
Banque8='IxodeDB'
Banque9='GlossinaDB'
Banque10='AutresDB'
Banque11='PrimatesDB'


id0=$(cut -c -3 ${File2}Banque0.txt)
id1=$(cut -c -3 ${File2}Banque1.txt)
id2=$(cut -c -3 ${File2}Banque2.txt)
id3=$(cut -c -3 ${File2}Banque3.txt)
id4=$(cut -c -3 ${File2}Banque4.txt)
id5=$(cut -c -3 ${File2}Banque5.txt)
id6=$(cut -c -3 ${File2}Banque6.txt)
id7=$(cut -c -3 ${File2}Banque7.txt)
id8=$(cut -c -3 ${File2}Banque8.txt)
id9=$(cut -c -3 ${File2}Banque9.txt)
id10=$(cut -c -3 ${File2}Banque10.txt)
id11=$(cut -c -3 ${File2}Banque11.txt)

	for Nomm in $Liste_noms
		do	
		dt=$('date')
		echo -e "	Debut soustraction ${Nomm}: $dt \n" >>log.txt

		#IDNom=$(echo "${ll}" | sed -e 's/.fasta//g')
		#echo "$IDNom"

		Nom=$(echo ${Nomm} | cut -c -2 )
		echo $Nom
	
		# *************** 
		#----------------Soustraction par alignement SR k4---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nomm}.fastq sur ${Banque0} \n"
			echo -e "-----------------------------------------------------------\n"

			echo ${Banque0} 

			# D1_no_Ndup.fastq
			`bwa mem -a -k 14 ${pathBanque}${Banque0}.fasta ${File}${Nomm}.fastq > ${Nom}${id0}-aln-se.sam` # different paramètres s
			#`bwa mem -a -k 14 ${Banque0}.fasta ${Nom}_my_Converted_D1_Test1_Int_k21.fastq > ${Nom}${id0}-aln-se.sam` # Paramètre utilisé  pour la db sous read/contig a k=14
			
			# -a  single-end 
			# -k = 19 minimum seed length  # longueur de k-mer  si diminu -> moins spécifique si augmente -> plus spécifique
			# -A = 1  Matching score : score des sequence matché
			# -c =10000 : nombre d'occurrence 	
			# -B = 4  penalty for a mismatch # pénalise les mismatch 
			# -O = 6 gap open penalty : ouverture # pénaliser les gaps
			# -E =1 gap extension penalty; a gap of size k cost {-O} + {-E}*k [1] # penalise les grand gaps
			# -L = 5  penalty for clipping : pénalité de coupure  # 			
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			#`samtools view -bS ${Nom}${id0}-aln-se.sam > ${Nom}${id0}.bam`	
			`samtools view -bS ${Nom}${id0}-aln-se.sam > ${Nom}${id0}.bam`		
			`samtools sort ${Nom}${id0}.bam ${Nom}${id0}_sorted`		
			`samtools index ${Nom}${id0}_sorted.bam`		
			`mv -f ${Nom}${id0}_sorted.bam.bai ${Nom}${id0}_sorted.bai`	
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}_sorted.bam > ${Nom}${id0}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}_sorted.bam > ${Nom}${id0}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}_sorted.bam > ${Nom}${id0}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}_sorted.bam > ${Nom}${id0}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}_Unmapped.bam > ${Nom}${id0}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}_Mapped.bam > ${Nom}${id0}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque0}.fasta ${Nom}${id0}_Mapped.bam | bcftools call -c -v >${Nom}${id0}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s ${Nom}${id0}Un.fastq ${Nom}${id0}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id0}_Resultat_Soustraction`
			rm *.sam 
			
			`mv ${Nom}${id0}_Unmapped.bam ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_Mapped.bam ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_Unmapped_stat.txt ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_Mapped_stat.txt ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_Mapped.vcf ${Nom}_${id0}_Resultat_Soustraction`
		
			rm *.bam 
			`cp ${Nom}${id0}Un.fastq ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}-aln-se.sam ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}.bam ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_sorted.bam ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_sorted.bai ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_sorted.sam ${Nom}_${id0}_Resultat_Soustraction`
			`mv ${Nom}${id0}_sorted_stat.txt ${Nom}_${id0}_Resultat_Soustraction`
	
	#SUPRIMER LES FICHIER.SAM 
	rm *.sam 


		#==============================================================================		
		
# samtools tview -d H G6_Wol_sorted.bam WolbachiaDB.fasta > R.html # Permet de visualiser l'alignement 
		
		#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un.fastq sur ${Banque1} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque1} #>Banque0.txt

			#`cp ${Nom}${id0}Un.fastq ${Nom}_Resultat_Soustraction`
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque1}.fasta ${Nom}${id0}Un.fastq > ${Nom}${id0}Un${id1}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}-aln-se.sam > ${Nom}${id0}Un${id1}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}.bam ${Nom}${id0}Un${id1}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}_sorted.bam.bai ${Nom}${id0}Un${id1}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}_sorted.bam > ${Nom}${id0}Un${id1}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}_sorted.bam > ${Nom}${id0}Un${id1}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}_sorted.bam > ${Nom}${id0}Un${id1}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}_sorted.bam > ${Nom}${id0}Un${id1}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}_Unmapped.bam > ${Nom}${id0}Un${id1}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}_Mapped.bam > ${Nom}${id0}Un${id1}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque1}.fasta ${Nom}${id0}Un${id1}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un.fastq  ${Nom}${id0}Un${id1}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id1}_Resultat_Soustraction`
			rm *.sam 
	
			`mv ${Nom}${id0}Un${id1}_Unmapped.bam ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_Mapped.bam ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_Unmapped_stat.txt ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_Mapped_stat.txt ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_Mapped.vcf ${Nom}_${id1}_Resultat_Soustraction`
			
			rm *.bam 
			`cp ${Nom}${id0}Un${id1}Un.fastq ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}-aln-se.sam ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}.bam ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_sorted.bam ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_sorted.bai ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_sorted.sam ${Nom}_${id1}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}_sorted_stat.txt ${Nom}_${id1}_Resultat_Soustraction`
		#==============================================================================
		
#**************************************************************************************************************************
			
		#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un.fastq sur ${Banque2} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque2} #>Banque0.txt

			#`cp ${Nom}${id0}Un.fastq ${Nom}_Resultat_Soustraction`
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque2}.fasta ${Nom}${id0}Un${id1}Un.fastq > ${Nom}${id0}Un${id1}Un${id2}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}.bam ${Nom}${id0}Un${id1}Un${id2}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque2}.fasta ${Nom}${id0}Un${id1}Un${id2}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}.fastq  ${Nom}${id0}Un${id1}Un${id2}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id2}_Resultat_Soustraction`
			rm *.sam 
			`mv ${Nom}${id0}Un${id1}Un${id2}_Unmapped.bam ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_Mapped.bam ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_Unmapped_stat.txt ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_Mapped_stat.txt ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_Mapped.vcf ${Nom}_${id2}_Resultat_Soustraction`
			
			rm *.bam 
			`cp ${Nom}${id0}Un${id1}Un${id2}.fastq ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}-aln-se.sam ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}.bam ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_sorted.bam ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_sorted.bai ${Nom}_${id2}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}_sorted.sam ${Nom}_${id2}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}_sorted_stat.txt ${Nom}_${id2}_Resultat_Soustraction`
	
	
		#==============================================================================
		# ***************






		
#**************************************************************************************************************************
			
		#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}.fastq sur ${Banque3} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque3} #>Banque0.txt

			#`cp ${Nom}${id0}Un.fastq ${Nom}_Resultat_Soustraction`
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque3}.fasta ${Nom}${id0}Un${id1}Un${id2}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}-aln-se.sam` # different paramètres 		

			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque3}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id3}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}.fastq ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted_stat.txt ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_Unmapped.bam ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped.bam ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_Unmapped_stat.txt ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped_stat.txt ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_Mapped.vcf ${Nom}_${id3}_Resultat_Soustraction`
			rm *.bam 
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}-aln-se.sam ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}.bam ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bam ${Nom}_${id3}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.bai ${Nom}_${id3}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}_sorted.sam ${Nom}_${id3}_Resultat_Soustraction`
			

		#==============================================================================
		
#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}.fastq sur ${Banque4} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque4} #>Banque0.txt

			#`cp ${Nom}${id0}Un.fastq ${Nom}_Resultat_Soustraction`
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque4}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque4}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id4}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.fastq ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}-aln-se.sam ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.bam ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bam ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.bai ${Nom}_${id4}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted.sam ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_sorted_stat.txt ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Unmapped.bam ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped.bam ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Unmapped_stat.txt ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped_stat.txt ${Nom}_${id4}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}_Mapped.vcf ${Nom}_${id4}_Resultat_Soustraction`
	
		#==============================================================================
		
#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.fastq sur ${Banque5} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque5} #>Banque0.txt

			#`cp ${Nom}${id0}Un.fastq ${Nom}_Resultat_Soustraction`
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque5}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque5}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id5}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.fastq ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}-aln-se.sam ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.bam ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bam ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.bai ${Nom}_${id5}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted.sam ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_sorted_stat.txt ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Unmapped.bam ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped.bam ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Unmapped_stat.txt ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped_stat.txt ${Nom}_${id5}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}_Mapped.vcf ${Nom}_${id5}_Resultat_Soustraction`
	
		#==============================================================================
	
#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.fastq sur ${Banque6} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque6} #>Banque0.txt

			#`cp ${Nom}${id0}Un.fastq ${Nom}_Resultat_Soustraction`
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque6}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque6}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id6}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.fastq ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}-aln-se.sam ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.bam ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bam ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.bai ${Nom}_${id6}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted.sam ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_sorted_stat.txt ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Unmapped.bam ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped.bam ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Unmapped_stat.txt ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped_stat.txt ${Nom}_${id6}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}_Mapped.vcf ${Nom}_${id6}_Resultat_Soustraction`
	
		#==============================================================================
		
#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.fastq sur ${Banque7} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque7} #>Banque0.txt

			
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque7}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque7}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id7}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.fastq ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}-aln-se.sam ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.bam ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bam ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.bai ${Nom}_${id7}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted.sam ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_sorted_stat.txt ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Unmapped.bam ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped.bam ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Unmapped_stat.txt ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped_stat.txt ${Nom}_${id7}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}_Mapped.vcf ${Nom}_${id7}_Resultat_Soustraction`
	
		#==============================================================================
		
#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.fastq sur ${Banque8} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque8} #>Banque0.txt
		
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque8}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque8}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id8}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.fastq ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}-aln-se.sam ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.bam ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bam ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.bai ${Nom}_${id8}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted.sam ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_sorted_stat.txt ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Unmapped.bam ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped.bam ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Unmapped_stat.txt ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped_stat.txt ${Nom}_${id8}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}_Mapped.vcf ${Nom}_${id8}_Resultat_Soustraction`
	
		#==============================================================================
		
#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.fastq sur ${Banque9} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque9} #>Banque0.txt
		
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque9}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque9}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id9}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.fastq ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}-aln-se.sam ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.bam ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bam ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.bai ${Nom}_${id9}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted.sam ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_sorted_stat.txt ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Unmapped.bam ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped.bam ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Unmapped_stat.txt ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped_stat.txt ${Nom}_${id9}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}_Mapped.vcf ${Nom}_${id9}_Resultat_Soustraction`
	
		#==============================================================================
		

#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.fastq sur ${Banque10} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque10} #>Banque0.txt
		
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque10}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque10}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id10}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.fastq ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}-aln-se.sam ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.bam ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bam ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.bai ${Nom}_${id10}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted.sam ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_sorted_stat.txt ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Unmapped.bam ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped.bam ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Unmapped_stat.txt ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped_stat.txt ${Nom}_${id10}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}_Mapped.vcf ${Nom}_${id10}_Resultat_Soustraction`
	
		#==============================================================================	

#----------------Soustraction par alignement---------------------------#
			echo -e "-----------------------------------------------------------\n"
			echo -e "Alignement de ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.fastq sur ${Banque11} \n"
			echo -e "-----------------------------------------------------------\n"
			echo ${Banque11} #>Banque0.txt
		
			# je fais varier la valeur du seed (k-mer) je la diminue pour êtres plus moins stringente et aligner plus de reads 
			`bwa mem -a -k 14 ${pathBanque}${Banque11}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}.fastq > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}-aln-se.sam` # different paramètres 		
			./Purge.sh
		#----------------Conversion sam en bam----------------------------------#
			`samtools view -bS ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}-aln-se.sam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}.bam`		
			`samtools sort ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}.bam ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted`		
			`samtools index ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam`		
			`mv -f ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam.bai ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bai`
			
		#---------------Trie des Mapped/Unmapped-------------------------------#
			#`samtools view -h ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.sam`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted_stat.txt`
			./Purge.sh
			`samtools view -b -f 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Unmapped.bam`
			./Purge.sh
			`samtools view -b -F 4 ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped.bam`
			./Purge.sh
		#---------------Statistiques des Mapped/Unmapped ----------------------#
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Unmapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Unmapped_stat.txt`
			./Purge.sh
			`samtools flagstat ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped.bam > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped_stat.txt`
			./Purge.sh
			`samtools mpileup -u -f ${pathBanque}${Banque11}.fasta ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped.bam | bcftools call -c -v > ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped.vcf`
			./Purge.sh
		#-------------Convertion .bam en .fastq----------------------#
			`samtools bam2fq -nO -s  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}.fastq  ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Unmapped.bam`
			./Purge.sh
		#-------Rangement les fichier resultats dans les repertoires correspondant------------#
			`mkdir ${Nom}_${id11}_Resultat_Soustraction`
			rm *.sam 
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}.fastq ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}-aln-se.sam ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}.bam ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bam ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.bai ${Nom}_${id11}_Resultat_Soustraction`
			#`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted.sam ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_sorted_stat.txt ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Unmapped.bam ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped.bam ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Unmapped_stat.txt ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped_stat.txt ${Nom}_${id11}_Resultat_Soustraction`
			`mv ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}_Mapped.vcf ${Nom}_${id11}_Resultat_Soustraction`
	
			`mkdir ${Nomm}_Resultat_Ass`
			`cp ${Nom}${id0}Un${id1}Un${id2}${id3}${id4}${id5}${id6}${id7}${id8}${id9}${id10}${id11}.fastq ${Nomm}_Resultat_Ass`
			`mv ${Nomm}_Resultat_Ass/ ${pathAssemblage}files`		

		#==============================================================================
		
		dt=$('date')
		echo -e "	Fin soustraction ${Nomm}: $dt \n" >>log.txt

		'date'		
		
		`mkdir ${Nomm}_Resultat`
		`mv *Resultat_Soustraction ${Nomm}_Resultat`
		#`mv ${Nomm}_Resultat/ /Volumes/Seagate\ Exp/Vanessa/BWA_test_param`
		
		#mv ${Nom}_Resultat_Test1/ /Volumes/Seagate\ Exp/Vanessa/BWA_test_param
		#mv D1_Aed_Resultat_Soustraction/ /Volumes/Seagate\ Exp/Vanessa/BWA_test_param
	done
fi 

dt=$('date')
echo -e "	Fin All soustraction : $dt \n" >>log.txt

'date'
