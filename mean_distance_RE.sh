#!/bin/bash

##### The following script calculates the average distance between ddRAD restriction enzyme cut-sites for each genome in a list #####
##### Brenna Levine - 28 October 2020 ####

#write help function to be printed if command line arguments are missing
helpFunction()
{
	echo ""
	echo -e "Usage: $0 -i <input_file> -r <P1_cutsite_seq> -f <P2_cutsite_seq>\n"
	echo -e "\t-i Input file containing genome file names and sequence lengths"
	echo -e "\t-r P1 cutsite recognition sequence (forward)"
	echo -e "\t-f P2 cutsite recognition sequence (forward)\n"
	echo -e "Contact Dr. Brenna Levine (levine.brenna.a@gmail.com) for help.\n"	    
	exit 1 #exit script after printing help
}

#use getopts to define command line arguments
while getopts "i:r:f:" opt;
do
	case ${opt} in
		i ) input_file="$OPTARG" ;;
		r ) P1_cutsite="$OPTARG" ;;
	 	f ) P2_cutsite="$OPTARG" ;;
		? ) helpFunction ;; #print helpFunction in case parameters are not provided
   	esac
done

#print helpFunction if any parameters are empty
if [ -z "$input_file" ] || [ -z "$P1_cutsite" ] || [ -z "$P2_cutsite" ]
then
	echo -e "\nEeeeek! Some or all required parameters are empty.\n";
	helpFunction
fi

#################################### Begin analysis #####################################

#echo status to STDOUT
echo -e "\n\n3...2..1...blast off!\n"

#echo parameter values
echo -e "\nThe following command line arguments have been passed as variables to the script:\n"
echo -e "\tinput file: '$input_file'"
echo -e "\tP1 cutsite: '$P1_cutsite'"
echo -e "\tP2 cutsite: '$P2_cutsite'"

################# Use tr to generate reverse cut-site recognition sequences

P1_cutsite_rev=$(echo $P1_cutsite | tr GCTA CGAT)
P2_cutsite_rev=$(echo $P2_cutsite | tr GCTA CGAT)

################# Echo forward and reverse cut site sequences
echo -e "\n\nP1 - Forward and reverse cut-site recognition sequences: $P1_cutsite $P1_cutsite_rev"
echo -e "P2 - Forward and reverse cut-site recognition sequences: $P2_cutsite $P2_cutsite_rev\n"

################# Assign genome file names to array elements

#initialize counter
COUNTER=0

#cut the first field and remove header from input file and write to temp file
cut -f 1 $input_file | grep -v Genome > genome_temp

#while loop through temp file to create array
while read line   #while reading the lines of the temp file
do
	genome_list[$COUNTER]=$line #assign the index of the array equal to COUNTER value to the value of the line in temp file
	let COUNTER=COUNTER+1 #increase COUNTER by one
done < genome_temp	#feed temp file to while loop

#echo some values to STDOUT to confirm genome file list matches genome array
echo -e "\nContents of input file: '$input_file'.\n" 
cat -n $input_file | column -t
echo -e "\nDoes the value of the genome file name in each line match the value of the genome file name in the corresponding array element?"

############### While loop to test whether the genome file contents and array elements match
#re-set counter
COUNTER=0

#while loop to test whether value of genome file name in line matches value of genome file name in array element
while read line #while reading the genome_temp file
do
	if [ $line = ${genome_list[$COUNTER]} ]; then   #if the value of line in genome_temp equals value of array element, then
		echo -e "\tYes."				#echo "Yes" to stdout
	else						#else
		echo -e "\tWARNING - Array element and line do not match!"	#echo warning to stdout
	fi						#close the if/else loop
	let COUNTER=COUNTER+1				#increase COUNTER variable by 1
done < genome_temp	#feed genome_temp file to while loop

#remove the temporary genome file
rm genome_temp

################# Assign genome file lengths to array elements

#re-set counter
COUNTER=0

#cut the second field and remove header from input file
cut -f 2 $input_file | grep -v Total  > genome_temp

#while loop through temp file to create array of genome lengths
while read line   #while reading the lines of the temp genome file
do
        genome_length[$COUNTER]=$line #assign the index of the array equalt to COUNTER value to the value of the line in temp genome file
        let COUNTER=COUNTER+1 #increase COUNTER variable by one
done < genome_temp      #feed genome_temp to while loop

#echo phrase
echo -e "\nDoes the value of the genome length in each line match the value of the genome length in the corresponding array element?"

#################### While loop to test whether the genome file contents and array match
#re-set counter
COUNTER=0

#while loop to test whether value of genome length in line matches value of genome length in array element
while read line #while reading line of temp file
do
        if [ $line = ${genome_length[$COUNTER]} ]; then 	#if the value of line in genome_temp equals value of array element, then
                echo -e "\tYes."     				#echo "Yes" to STDOUT
        else                    				#else
                echo -e "\tWARNING - Array element and line do not match!"      #echo warning to STDOUT
        fi                      				#close the if/else loop
        let COUNTER=COUNTER+1   				#increase the counter by 1
done < genome_temp      #feed genome_temp file to while loop

#remove temporary file
rm genome_temp

#################### While loop through arrays to calculate P1/P2 cut-site frequencies
#re-set counter
COUNTER=0

#calculate number of elements in genome_list and genome_length arrays
NUM_LIST=${#genome_list[@]}
NUM_LENGTH=${#genome_length[@]}

#confirm that genome_list and genome_length arrays have the same number of elements
if [ $NUM_LIST = $NUM_LENGTH ]; then 	#if the two arrays have the same number of elements, then
	echo -e "\n\nOK - The genome_list and genome_length arrays have the same number of elements: $NUM_LIST"
else	
	echo -e "\n\nWARNING: The genome_list and genome_length arrays do not have the same number of elements!" 
fi

#create headers in two new files
echo "P1_Distance" > P1_freq.txt	#echo header to temporary file
echo "P2_Distance" > P2_freq.txt	#echo header to temporary file

#echo status to STDOUT
echo -e "\n\nCalculating cut-site frequencies for genome files in list............\n"

while [ $COUNTER -lt $NUM_LIST ]	#while the value of COUNTER is less than the value of NUM_LIST
do
	echo -e "..........................................................\n" #echo line of dots
	num_P1_for=$(grep -oi $P1_cutsite ${genome_list[$COUNTER]} | wc -l) #number of P1 forward cutsites = grep sequence from genome file and run line count
	echo -e "..........................................................\n" #echo line of dots
	num_P1_rev=$(grep -oi $P1_cutsite_rev ${genome_list[$COUNTER]} | wc -l)  #number of P1 reverse cutsites = grep sequence from genome file and run line count
	echo -e "Number of forward P1 cut-sites in ${genome_list[$COUNTER]}: $num_P1_for\n"
	echo -e "Number of reverse P1 cut-sites in ${genome_list[$COUNTER]}: $num_P1_rev\n"
	P1_total=$( expr $num_P1_for + $num_P1_rev) #add forward and reverse totals to get total number of P1 cut-sites in genome
	echo -e "Total number of P1 cut-sites in ${genome_list[$COUNTER]}: $P1_total\n"
	P1_freq=$( expr ${genome_length[$COUNTER]} / $P1_total)

	echo -e "..........................................................\n" #echo line of dots
        num_P2_for=$(grep -oi $P2_cutsite ${genome_list[$COUNTER]} | wc -l) #number of P2 forward cutsites = grep sequence from genome file$
        echo -e "..........................................................\n" #echo line of dots
        num_P2_rev=$(grep -oi $P2_cutsite_rev ${genome_list[$COUNTER]} | wc -l)  #number of P2 reverse cutsites = grep sequence from genome$
        echo -e "Number of forward P2 cut-sites in ${genome_list[$COUNTER]}: $num_P2_for\n"
        echo -e "Number of reverse P2 cut-sites in ${genome_list[$COUNTER]}: $num_P2_rev\n"
        P2_total=$( expr $num_P2_for + $num_P2_rev) #add forward and reverse totals to get total number of P2 cut-sites in genome
        echo -e "Total number of P2 cut-sites in ${genome_list[$COUNTER]}: $P2_total\n"
        P2_freq=$( expr ${genome_length[$COUNTER]} / $P2_total)

	echo -e "Mean number of base pairs between P1 cut-sites in ${genome_list[$COUNTER]}: $P1_freq\n"
	echo -e "Mean number of base pairs between P2 cut-sites in ${genome_list[$COUNTER]}: $P2_freq\n"
	echo -e "Writing totals to temporary files..........................\n"

	#write totals to two temporary files
	echo $P1_freq >> P1_freq.txt	#write total number of P1 fragments to a temporary file
	echo $P2_freq >> P2_freq.txt	#write total number of P2 fragments to a temporary file

	let COUNTER=COUNTER+1	#increase COUNTER variable by 1
done

#paste original genome list file to cut-site frequency files to produce final output file
paste $input_file P1_freq.txt P2_freq.txt > mean_distance_out.txt

#remove intermediate files
rm P1_freq.txt P2_freq.txt

#echo final phrase to STDOUT
echo -e "\n\nFor each genome file, mean distances between cut sites are reported in mean_distance_out.txt.\n"
echo -e "\nSee you next time!\n\n\n"
