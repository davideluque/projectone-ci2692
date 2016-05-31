# Universidad Simón Bolívar
# CI2692: Laboratory of Algorithms and Structures II

# Project 1: Implementing of algorithms to emulate the synthesis process of 
# 			 proteins.
# Authors: David Cabeza 13-10191, Fabiola Martinez 13-10838
# Last edit: Sun 22, May 2016
import sys
from colorama import *
init()

# Function get_complemet: get a simple strand and return its complement
# Parameters: 
# 	Sequence: a simple strand of DNA
def get_complement(sequence):
	complement = ""
	for nucleobase in sequence:
		if nucleobase == "A":
			complement += "T"
		elif nucleobase == "T" or nucleobase == 'U':
			complement += "A"
		elif nucleobase == "G":
			complement += "C"
		elif nucleobase == "C":
			complement += "G"

	return complement

# Function read strands: read a file that it conteins strands of DNA, then the creates the complement and finally write doble DNA in pairs
# Parameters:
#
def read_strands(DnaFile):
	DNAsimple = []
	with open(DnaFile, 'r') as f:
		for line in f:
			line = line.rstrip()
			DNAsimple.append(line)
		print('Beginning to read a strand file...')
		print(DNAsimple)	
	f.close()

	complements = []
	for elements in DNAsimple:
		get_complement(elements)
		complements.append(get_complement(elements))
	print('Getting complement...')
	print(complements)

	DNADoble = []
	for i in range(0,len(complements)):
		for core in complements[i]:
			if core == "A":
				DNADoble += "AT"
			elif core == "T":
				DNADoble += "TA"
			elif core == "G":
				DNADoble += "GC"
			elif core == "C":
				DNADoble += "CG"
	print('Getting DNADoble...')
	print(DNADoble)
	
	print('Create the pars...')
	for i in range(0, len(DNADoble), 2):
		print('(' + Fore.RED + DNADoble[i], Fore.BLUE + DNADoble[i+1], Style.RESET_ALL+')', sep ='', end='')

print(read_strands('data.txt'))


class ADNDoble(object):
	"""Class whose instances are sequences of simple ADN. This class manipulates 
	its associated double sequence"""
	def __init__(self, simpleseq):
		self.simpleseq = simpleseq.upper()
		self.sseqlenght = len(self.simpleseq)
		self.complement = ""
		self.doubleseq = ""

# Function zip: get a simple strand, then mades the complement and finally creates the doble DNA
# Parameters: 
# 
	def zip(self):
		print('\nBeginning to zip the sequence', ADN.simpleseq)
		for nucleobase in self.simpleseq:
			if nucleobase == "A":
				self.doubleseq += "AT"
			elif nucleobase == "T":
				self.doubleseq += "TA"
			elif nucleobase == "G":
				self.doubleseq += "GC"
			elif nucleobase == "C":
				self.doubleseq += "CG"

		self.complement = get_complement(self.simpleseq)
		print('The complement for the DNA sequence given is:', self.complement)
		print('The double sequence for the given is:', self.doubleseq)
		print('Zip process was made successfully!!')
	
# Function unzip: get a doble strand and divide it in two simple strand
# Parameters: 
# 	
	def unzip(self):
		print('\nBeginning to unzip the sequence', ADN.doubleseq)	
		seqOne = ""
		for i in range(0, len(self.doubleseq), 2):
			seqOne += self.doubleseq[i]
		seqTwo = ""
		for j in range(1, len(self.doubleseq), 2):
			seqTwo += self.doubleseq[j]
		print('The simple two sequences are:')
		print(seqOne,'and', seqTwo)

# Function mitosis: get a doble strand then duplicate this and return two doble stra
# Parameters: 
# 
	def mitosis(self):
		print('\nBeginning the mitosis process for the sequence', ADN.doubleseq)	
		print('Beginning to unzip the double sequence')
		seqOne = ""
		for i in range(0, len(self.doubleseq), 2):
			seqOne += self.doubleseq[i]
		seqTwo = ""
		for j in range(1, len(self.doubleseq), 2):
			seqTwo += self.doubleseq[j]
		print('Beginning to calculate the complements')
		CseqOne = get_complement(seqOne)
		CseqTwo = get_complement(seqTwo)

# Function buscar: get a simple substrand and search its match in doble strand
# Parameters: 
# 

	def buscar(self, subcadena):
		print("\nBeginning search subsequence process..")
		complement = get_complement(subcadena)
		if subcadena in self.doubleseq:
			print('MATCH, we found', subcadena, 'in', self.doubleseq)
		if complement in self.doubleseq:
			print('MATCH, we found the complement', complement, 'of', subcadena,'in', self.doubleseq)

# Function imprimir: write the strand in par form
# Parameters: 
#
	def imprimir(self):
		print()
		for i in range(0, len(self.doubleseq), 2):
			print('(' + Fore.RED + self.doubleseq[i], Fore.BLUE + self.doubleseq[i+1], Style.RESET_ALL+')', sep ='', end='')
		print()

# Function write: it allows write in a file and save the proteins or the strands
# Parameters:
#

	def write(self, file):
		if not '.txt' in file:
			file += '.txt'
		with open(file, 'a') as f:
			f.write(self.doubleseq)

class ADNSimple(object):
	"""Instances: Simple ADN sequence"""
	def __init__(self, sequence):
		self.sequence = sequence
		self.complim = ""

# Function complimet: get a simple strand and return its complement
# Parameters: 
#

	def complement(self):
		for nucleobase in self.sequence:
			if nucleobase == "A":
				self.complim += "T"
			elif nucleobase == "T":
				self.complim += "A"
			elif nucleobase == "G":
				self.complim += "C"
			elif nucleobase == "C":
				self.complim += "G"

# Function tansilerate: change the Timina by Uracilo
# Parameters: 
#

	def transliterate(self):
		if not self.complim:
			self.complim == self.complement()

		self.ARNt = ""
		for nb in self.complim:
			if nb == "T":
				self.ARNt += "U"
			else:
				self.ARNt += nb

		return self.ARNt

class ARNt(ADNSimple):
	"""docstring for ARNt"""
	def __init__(self, sequence):
		self.sequence = sequence
		self.size = len(sequence)
		self.aminoacids = ['Phe', 'Phe', 'Leu','Leu', 'Ser', 'Ser', 'Ser', 
		'Ser', 'Tyr', 'Tyr', 'STOP', 'STOP', 'Cys', 'Cys', 'STOP', 'Trp',
		'Leu', 'Leu','Leu','Leu', 'Pro', 'Pro', 'Pro', 'Pro', 'His', 'His',
		'Gln', 'Gln', 'Arg', 'Arg', 'Arg', 'Arg', 'Ile', 'Ile', 'Ile', 'Met',
		'Thr', 'Thr', 'Thr', 'Thr', 'Asn', 'Asn', 'Lys', 'Lys', 'Ser', 'Ser',
		'Arg', 'Arg', 'Val', 'Val', 'Val', 'Val', 'Ala', 'Ala', 'Ala', 'Ala', 
		'Asp', 'Asp', 'Glu', 'Glu', 'Gly', 'Gly', 'Gly', 'Gly']
		self.proteins = []
		self.letras = ['U', 'C', 'A', 'G']
		self.ARNJunks = []

# Function translate: search the start codon and process the secuence until find the stop protein 
# Parameters:
	def translate(self):
		start = False
		for i in range(0, self.size, 3):
			trio = self.sequence[i]+self.sequence[i+1]+self.sequence[i+2]
			if i == 69:
				print(trio)
				print(start)
			letras = False
			if not start and trio == 'AUG':
				start = True
				start_index = i
				end = False
				letras = True
				temp = []
			if start and trio != 'AUG' and trio != 'UAA' and trio != 'UAG'\
			and trio != 'UGA':
				value = self.letras.index(self.sequence[i]) * 16 \
				+ self.letras.index(self.sequence[i+1]) * 4 \
				+ self.letras.index(self.sequence[i+2]) * 1
				temp.append(self.aminoacids[value])
				letras = True
			if not end and trio == 'UAA' or trio == 'UAG' or trio == 'UGA':	
				self.proteins.append(temp)
				del temp
				start = False
				end = True
				letras = True
			if start and not end and not letras:
				j = start_index + 3
				temporary = []
				while j != i:
					junk = self.sequence[j] + self.sequence[j+1] + self.sequence[j+2]
					j += 3
					temporary.append(junk)
				self.ARNJunks.append(temporary)
				del temporary
				start = False
				end = True
				del temp
		
		print(self.proteins)
		print(self.ARNJunks)

# Main program			
if __name__ == '__main__':
	ADN = ADNDoble("ATGGGCAATCGGTTTGC")
	ADN.zip()
	ADN.unzip()
	ADN.buscar("ATA")
	ADN.mitosis()
	ADN.imprimir()

	var = get_complement('ATGTTTTTCTTATTGTCTTCCTCATCGTATTACTAAATGACGATAGTAGATTGAATGTTCTAAATGTTTATGTCTTAA')
	ADNs = ADNSimple(var)
	ADNs.complement() # If you don't call this method first, the transliterate method
	transliterated = ADNs.transliterate() # will do it automatically.

	print(transliterated)
	protein = ARNt(transliterated)
	protein.translate()
	protein.size()