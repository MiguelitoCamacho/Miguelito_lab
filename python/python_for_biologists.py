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