# This is the training course Frank found on internet. https://www.pythonforbiologists.org/

# Introduction to Python 

# 1. Welcome. 

print ( "Welcome to Python for biologists!")

# 2. Beginnins 
# References:
# https://docs.python.org/3/tutorial/introduction.html#an-informal-introduction-to-python
# https://www.ncbi.nlm.nih.gov/nuccore/U00096
# https://en.wikipedia.org/wiki/Chargaff%27s_rules

organism = "Escherichia coli"
strain = ' stre. K-12 substr. MG1665'
print ("DEFINITION: " + organism + " " + strain) 

print (organism)

number_of_bsp = 4641642
print( "Number of base paris:", number_of_bsp)

percent_A = 24.7 
percent_T = 23.6
percents_AGCT = [ percent_A, 26.0, 25.7, percent_T]
 
print( "(A, G, C, T) =", percents_AGCT)

ratio_AT = percent_A / percent_T 
ratio_GC = percents_AGCT[1]/ percents_AGCT[2]

print("ratio_ AT:", ratio_AT, "ratio G/C", ratio_GC)

E_coli = (organism, ratio_AT, ratio_GC)
print(E_coli)

# 3. print 

human_genes = 20000
print('You have', human_genes, 'genes')

print( 'You have', end = ' ')
print(human_genes, end = '?')
print('genes')

print( 'You have', 
      human_genes,
        "genes")

print( 'You have ' + '\'' + str(human_genes) + '\'' + '  genes')

import math

print(math.pi)

print( 'The value of pi is %10s %5.3f' %( '--->', math.pi))

print (chr(7))

# 4. File name: buggy.py
#
#This short Python code contains a number of interntional bugs. Correct
# them with the help of the error messages.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://americanhistory.si.edu/collections/search/object/nmah_334663
# https://en.wikipedia.org/wiki/Software_bug
# ------------------------------------------------------------------
 
#primero el ejemplo estaba como human_genes = 20,000. Then the instructor explained python does not like comma in numbers. 
human_genes = 20000
 
protein_name = "GFP" ;
 
print('You have', human_genes, 'genes')
print(protein_name, "stands for green fluorescent protein")


# exercises 
# 1. ASCII art: Write a Python program that prints out the image below.

print(" __")              
print("|  \\   |\\ |    /\\ ") 
print("|__/   | \\|   /--\\")
print("                ___             ___           ___    ____      ")
print( " |\\   /|   |   | ___   |   |   |___  |     |   |    |    | ")
print(" | \\ / |   |   |____|  |___|   |___  |___  |   |    |____|     ")

# 2. Finding bugs: Find and exterminate the bugs in the Python code below

# Please correct my errors. 
first_10_bp = "gggtgcgacg"
second_10_bp = "attcattgtt"
gene = first_10_bp + second_10_bp

print('The first 20 bp of BRAC2 gene are', gene)

#so that it prints out:
#The first 20 bp of BRAC2 gene are gggtgcgacgattcattgtt

# 3.When precise control with the output of your program is needed (lining up columns, decimal points, etc.),
#    the print function of Python is your tool. It takes two arguments:


print("%15s %10s %10s" % ("Amino acid", "1-letter", "codon"))
print("%15s %10s %10s" % ("----------", "--------", "-----"))
print("%15s %10s %10s" % ("Serine", "S", "AGT"))
print("%15s %10s %10s" % ("", "", "AGC"))
print("%15s %10s %10s" % ("Cysteine", "C", "TGT"))

  
# Data types 

#Numeric 
import math
human_genes = 21306
US_population = 328918373
print ('Number of human genes:', human_genes)
print ('Number of human genes in US:', US_population)

exon_per_gene = 8.9 
print('Human exons per gene:', exon_per_gene)

human_exons = human_genes * exon_per_gene 

print ('Humber of human exons:', human_exons)

human_exons = int (human_exons)

print ("Approximate number of human exon:", human_exons)

firstProduct = (9.4 * 0.2321) * 5.6
secondProduct = 9.4 * (0.2321 * 5.6)

print (firstProduct - secondProduct)

two_pie = 2.0 * math.pi
print (' Two_pi = ', two_pie)
print ("sin(two_pie) =", math.sin(two_pie))


# Strings 

protein = 'GFP'
protein_seq_begin = 'MSKGEELFTG'
protein_seq_end = 'HGMDELYK'
protein_seq = protein_seq_begin + '...' + protein_seq_end
print( protein_seq)

print ( 'Protein sequence of GPF: ' + protein_seq)
print ( 'Protein sequence of GPF:',protein_seq)

DNA_seq = 'atgagtaaag...actatacaaa'
DNA_seq = DNA_seq.upper()
print ('DNA Sequence: ' + DNA_seq)

print ( 'The second nucleotide: ' + DNA_seq[1])
print ( 'The last nucleotide: ' + DNA_seq[-1])

first_codon = DNA_seq[0:3]  # it could be also [ : 3]
last_codon = DNA_seq[-3:]
print( 'The first codon: ' + first_codon)
print ('The last codon: ' + last_codon)

# Lists 

stop_codons = [ 'TAA', 'tAG']
print( stop_codons)

first_stop_codon = stop_codons[0]
print(first_stop_codon)

stop_codons [1] = 'TAG'
print( stop_codons)

stop_codons.append('TGA') 
print( stop_codons)

number_of_stop_codons = len (stop_codons)
print ('There are', number_of_stop_codons, 'stop codons')

DNA_seq = ''.join(stop_codons)
print( DNA_seq)

DNA_list = list(DNA_seq)
print(DNA_list)
                              

second_codon = DNA_list[3:6]
print('This is the second codon: ', second_codon)

print (''.join(second_codon))


DNA_list_duplicate = DNA_list.copy()
print (DNA_list_duplicate)

DNA_list_duplicate.insert(5, '?')
print (DNA_list_duplicate)


DNA_list_duplicate.pop(5)
print(DNA_list_duplicate)

# Tupples 

Histidine = ('H', 'CAT', 'CAC')
Lysine = 'K', 'AAA', 'AAG'
Arginine = ( 'R', 'CGT', 'CGA', 'CGG', 'AGA', 'AGG')

print ('Histidine', Histidine)
print ('Lysine', Lysine)
print('Arginine', Arginine)

basic = [ Histidine, Lysine]
print ('Basic amino acids:', basic)
basic.append(Arginine)
print ('Basic amino acids:', basic)

His = basic[0]
print('His:', His)

His_codons = basic[0][1:]

print('His codons:', His_codons)

codon1, codon2 = His_codons
print('codon1:', codon1)
print('codon2:', codon2)


protein_seq = basic[0][0] + basic[1][0] + basic[2][0]

print('Protein seq:', protein_seq)

#Dictionary 

restriction_enzymes = { 'EcoRI': 'GAATTC',
                        'AluI': 'AGCT',
                        'NotI' : 'GCGGCCGC',
                        'TaqI': 'TCGA'
  }
print (restriction_enzymes)

# how we extract keys and values and make them list 

keys = list( restriction_enzymes.keys())

print( 'Keys as a list:', keys)

values = list(restriction_enzymes.values())
print('Values as a list:', values)

mykey = 'crispr'
check = mykey in restriction_enzymes
print (' Is', mykey, 'key in the dcitionary?', check)

EcorRI_value = restriction_enzymes['EcoRI']
print ('The recognition site of EcoRI is', EcorRI_value)

restriction_enzymes['EcoRV'] = 'GATATC'
print('With the new item:', restriction_enzymes)

del restriction_enzymes['EcoRV']
print('Original dictionary:', restriction_enzymes)

# Exercises II 
#1. Arithmetic puzzle: Create two numeric float variables
x = 0.1 * 10.0
y = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1
print(x-y)

#2. BRAC2 mRNA: Visit the Web Page. 

BRAC2_DNA_seq = 'gggtgcgacg attcattgtt ttcggacaag tggataggca accactaccg gtggattgtc'
print(BRAC2_DNA_seq)
BRAC2_DNA_seq = BRAC2_DNA_seq.replace(' ','')
print(BRAC2_DNA_seq)
Upper_case_BRAC2 = BRAC2_DNA_seq.upper()
lower_case_BRAC2 = BRAC2_DNA_seq.lower()
print('This is the BRAC2 se in upper case', Upper_case_BRAC2)
print('This is the BRAC2 seq in lowercase', lower_case_BRAC2)

#3. List of codon acidi acids: Creat a list with the DNA codons for Glutamine

codons = [ 'CAA', 'CAG']

codons.extend(['GAU', 'GAC'])

print(codons)


print("Glutamine codons:", codons[:2])
print('Aspartic acid cordons:', codons[3:])


# 4. append() and pop(): Sometimes it is convenient to add to or subtract from a list without worrying about indecies.

codons.append( 'GAT')
last_element = codons.pop()

print( last_element)

#5. Sorting a list of string 

sorted_list = sorted(codons)
print (' This is the sorted list', sorted_list)

numbers = [21, 13, 0, -3.1, 21, 33, 2]
numbers.sort()
print ('This are the sorted numbers:', numbers)

#7. Dictionary of codons for acidic acids 

import operator
codon_to_amino_acid = {
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine',
    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine',
    'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cysteine', 'UGC': 'Cysteine', 'UGA': 'Stop', 'UGG': 'Tryptophan',
    'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    'CAU': 'Histidine', 'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine',
    'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine',
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine', 'AUG': 'Methionine',
    'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    'AAU': 'Asparagine', 'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine',
    'AGU': 'Serine', 'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine',
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    'GAU': 'Aspartic acid', 'GAC': 'Aspartic acid', 'GAA': 'Glutamic acid', 'GAG': 'Glutamic acid',
    'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'
}


#9 Sorting a dictionanry by keys or by values: 

codon_to_amino_acid_sorted_by_kets = sorted(codon_to_amino_acid.items(), key =operator.itemgetter(0))
#print(codon_to_amino_acid_sorted_by_kets)

codon_to_amino_acid_sorted_by_values = sorted(codon_to_amino_acid.items(), key=operator.itemgetter(1))
#print(codon_to_amino_acid_sorted_by_values)

# Length 

zika_DNA = 'AGTTGTTGATCTGTGT'
zika_DNA_length = len(zika_DNA)
print('The first', zika_DNA_length, 'nucleotides', 'fo zika virus DNA are', zika_DNA)

# Contatenation 

GFP_seq = 'MSKGEELFTG...HGMDELYK '
print('Green flourescnece protein:', GFP_seq)

M_codon = 'AUG'
S_codon = 'UCA'
K_codon = 'AAA'
G_codon = 'GGU'

RNA_seq = M_codon 
RNA_seq = RNA_seq + S_codon
print('RNA sequence', RNA_seq)

RNA_seq = RNA_seq + K_codon + G_codon 
print(RNA_seq, 'could code amino acid sequence MSKG')

# Slicing

human = 'TTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTTACCCATCAACAACCGCTATGTATT'
print('Human D-loop:', human)

mycodon = 'CAT'

index_mycodon = human.find(mycodon)
print ('First', mycodon, 'index:', index_mycodon)

first_codon = human[index_mycodon + 3: index_mycodon + 6]

print('First codon after', mycodon, ':', first_codon)

second_codon = human[index_mycodon + 6: index_mycodon + 9]
print('Second codon after', mycodon, ':', second_codon)

next_to_last_codon = human[ -6: -3]
print('Next to last codon:', next_to_last_codon)

last_codon = human[-3:]
print('Last codon:', last_codon)



2 == (1 + 1)