#Part 1 script
#!/bin/bash

#print your name
echo 'Mariam Mohamed'
#create folder with your name
mkdir mariam
#make another directory titled 'biocomputing' and change to that directory with one line.
mkdir biocomputing && cd biocomputing
#download 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
#move the fna files to the directory with your name
mv *.fna ../mariam
#remove the duplicate gbk file
rm *.gbk.*
#confirm if fna is mutant or wild
cd ../

if grep -q "tatata" mariam/wildtype.fna; then
  echo "Mutant file"
  echo "Mutant file" >> results.txt
elif grep -q "tata" mariam/wildtype.fna; then
  echo "Wild type file"
  echo "Wild type file" >> results.txt
fi

#count lines without header
lines_with_header=$(wc -l biocomputing/wildtype.gbk | awk '{print $1}')
lines=$((lines_with_header -1))
echo "there are $lines in the gbk file"

#print the length of the sequence using 'LOCUS' 
grep 'LOCUS' biocomputing/wildtype.gbk | awk '{print $3}'
# Print the source organism of the .gbk file. (Use the SOURCE tag in the first line)
grep 'SOURCE' biocomputing/wildtype.gbk | awk '{print $2, $3}'

# List all the gene names of the .gbk file. Hint {grep '/gene='}

grep '/gene' wildtype.gbk 
# Clear your terminal space and print all commands used today
clear
#List the files in the two folders
ls mariam/ biocomputing/





