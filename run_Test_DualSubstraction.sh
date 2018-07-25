#!/bin/bash

if [ -d "result" ]
	then
		echo "result existe !"
		rm -R result
		mkdir result
	else
		echo " result no existe "
		mkdir result
fi

pth=$(pwd)/


echo -e "-------------------------------------------------------------------"
echo -e "   	 	0  - Index les banques de soustraction    "
echo -e "-------------------------------------------------------------------"

#bash ./Index.sh


echo -e "-------------------------------------------------------------------"
echo -e "          1  - Soustraction Reads       		 "
echo -e "-------------------------------------------------------------------"

bash ./Substraction_SR.sh


echo -e "-------------------------------------------------------------------"
echo -e "            2  - Assemblage       "
echo -e "-------------------------------------------------------------------"


bash ./Assembly.sh

echo -e "-------------------------------------------------------------------"
echo -e "            3  - Soustraction Contigs     "
echo -e "-------------------------------------------------------------------"


bash ./Substraction_SC.sh


`mv result_SR result_SC result_Assembly result`
`mv data_SRSC dataAssembly data_Ass_SR result`