# Universidad Simón Bolívar
# CI2692: Laboratory of Algorithms and Structures II
# Project 1: Implementing of algorithms to emulate the synthesis process of 
# 			 proteins.
# Authors: David Cabeza 13-10191, Fabiola Martinez 13-10838
# Last edit: Mon 30, May 2016

import sys
from random import randint

################### get_complement_tRNA METHOD #################################
# Description: This method receives a tRNA sequence and complements it changing
#			   Uracil for Adenine, Guanine for Cytosine and vice versa.
# Input: Receives a sequence of tRNA to be complemented 
# Output: Returns the complement of the given tRNA sequence
################################################################################
def get_complement_tRNA(sequence):
	complement = ""
	for base in sequence:
		if base == "A":
			complement += "U"
		elif base == "G":
			complement += "C"
		elif base == "U":
			complement += "A"
		else:
			complement += "G"

	return complement

######################## get_complement METHOD #################################
# Description: This method receives a DNA sequence and complements it changing
#			   Timine for Adenine, Guanine for Cytosine and vice versa.
# Input: Receives a sequence of DNA to be complemented 
# Output: Returns the complement of the given DNA sequence
################################################################################
def get_complement(sequence):
	complement = ""
	for nucleobase in sequence:
		if nucleobase == "A":
			complement += "T"
		elif nucleobase == "T":
			complement += "A"
		elif nucleobase == "G":
			complement += "C"
		elif nucleobase == "C":
			complement += "G"

	return complement

def Partition(A, p, r):
	x = len(A[r])
	i = p - 1
	# Invariant
	assert(all(len(A[i]) >= x for i in range(0, i+1)))
	for j in range(p, r):
		if len(A[j]) > x:
			i += 1
			A[i], A[j] = A[j], A[i]
		elif len(A[j]) == x:
			if min(A[j]) <= min(A[r]):
				i+= 1
				A[i], A[j] = A[j], A[i]
		# Invariant
		assert(all(len(A[i]) >= x for i in range(0, i+1)))
	A[i+1], A[r] = A[r], A[i+1]
	return i + 1

def quicksort_randomizedPartition(A, p, r):
	i = randint(p, r)
	A[i], A[r] = A[r], A[i]
	return Partition(A, p, r)

def quicksort(A, p = None, r = None):
	p = p if p else 0 
	r = r if r else len(A) - 1

	if p < r:
		q = quicksort_randomizedPartition(A, p, r)
		quicksort(A, p, q-1)
		quicksort(A, q+1, r)

	return A

def heapsort(A):
	n = len(A)	
	buildHeap(A)

	for i in range(n-1, 0, -1):
		A[0], A[i] = A[i], A[0]
		heapify(A, 0, i)

def buildHeap(A):
	n = len(A)
	for i in range((n//2), -1, -1):
		heapify(A, i, n)

def heapify(A, idx, mx):
	largest = idx
	left = (2*idx) + 1
	right = (2*idx) + 2

	if left < mx and A[left][1] < A[idx][1]:
		largest = left
	if right < mx and A[right][1] < A[largest][1]:
		largest = right
	if largest != idx:
		A[idx], A[largest] = A[largest], A[idx]
		heapify(A, largest, mx)

class DNADoble(object):
	"""Class whose instances are sequences of simple DNA. This class manipulates 
	its associated double sequence"""
	def __init__(self, simpleseq):
		self.simpleseq = simpleseq.upper()
		self.sseqlenght = len(self.simpleseq)
		self.complement = ""
		self.doubleseq = ""

	def zip(self):
		print('\nBeginning to zip the sequence', DNA.simpleseq)
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
	
	def unzip(self):
		print('\nBeginning to unzip the sequence', DNA.doubleseq)	
		seqOne = ""
		for i in range(0, len(self.doubleseq), 2):
			seqOne += self.doubleseq[i]
		seqTwo = ""
		for j in range(1, len(self.doubleseq), 2):
			seqTwo += self.doubleseq[j]
		print('The simple two sequences are:')
		print(seqOne,'and', seqTwo)

	def mitosis(self):
		print('\nBeginning the mitosis process for the sequence', DNA.doubleseq)	
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

	def buscar(self, subcadena):
		print("\nBeginning search subsequence process..")
		complement = get_complement(subcadena)
		if subcadena in self.doubleseq:
			print('MATCH, we found', subcadena, 'in', self.doubleseq)
		if complement in self.doubleseq:
			print('MATCH, we found the complement', complement, 'of', subcadena,'in', self.doubleseq)

	def imprimir(self):
		print()
		for i in range(0, len(self.doubleseq), 2):
			print('(',self.doubleseq[i],self.doubleseq[i+1],')', sep ='', end='')
		print()

	def write(self, file):
		if not '.txt' in file:
			file += '.txt'
		with open(file, 'a') as f:
			f.write(self.doubleseq)

class DNASimple(object):
	"""Instances: Simple DNA sequence"""
	def __init__(self, sequence):
		self.sequence = sequence
		self.complem = ""

	def complement(self):
		for nucleobase in self.sequence:
			if nucleobase == "A":
				self.complem += "T"
			elif nucleobase == "T":
				self.complem += "A"
			elif nucleobase == "G":
				self.complem += "C"
			elif nucleobase == "C":
				self.complem += "G"

	def transliterate(self):
		if not self.complem:
			self.complem == self.complement()

		self.tRNA = ""
		for nb in self.complem:
			if nb == "T":
				self.tRNA += "U"
			else:
				self.tRNA += nb

		return self.tRNA

class tRNA():
	"""docstring for tRNA"""
	def __init__(self, sequence):
		self.sequence = sequence
		self.size = len(self.sequence)
		self.complement = ""
		self.complementsize = len(self.complement)
		self.bases = ['U', 'C', 'A', 'G']
		self.aminoacids = ['Phe', 'Phe', 'Leu','Leu', 'Ser', 'Ser', 'Ser', 
		'Ser', 'Tyr', 'Tyr', 'STOP', 'STOP', 'Cys', 'Cys', 'STOP', 'Trp',
		'Leu', 'Leu','Leu','Leu', 'Pro', 'Pro', 'Pro', 'Pro', 'His', 'His',
		'Gln', 'Gln', 'Arg', 'Arg', 'Arg', 'Arg', 'Ile', 'Ile', 'Ile', 'Met',
		'Thr', 'Thr', 'Thr', 'Thr', 'Asn', 'Asn', 'Lys', 'Lys', 'Ser', 'Ser',
		'Arg', 'Arg', 'Val', 'Val', 'Val', 'Val', 'Ala', 'Ala', 'Ala', 'Ala', 
		'Asp', 'Asp', 'Glu', 'Glu', 'Gly', 'Gly', 'Gly', 'Gly']
		self.proteins = []
		self.tRNATrash = []
		self.DNATrash = []
		self.frequencies = []

	#################### Translate METHOD ####################################
	# Description: This method translates a sequence into proteins, checks if
	# exists 
	##########################################################################
	def Translate(self):
		self.complement = get_complement_tRNA(self.sequence)
		# Cambiar todos los self.sequence y self.size por los del complemento
		#print(self.sequence)
		Start = False
		End = True
		for i in range(0, self.size, 3):
			Trio = self.sequence[i] + self.sequence[i+1] + self.sequence[i+2]
			if not self.frequencies:
				self.frequencies.append([Trio, 1])
			else:
				j = 0
				while j < len(self.frequencies):
					if Trio == self.frequencies[j][0]: 
						self.frequencies[j][1] += 1
						break
					j += 1
				else:	
					self.frequencies.append([Trio, 1])
			if Trio == 'AUG' and not Start:
				Start = True
				Start_index = i
				End = False
				Temp_proteins = []
				continue
			if Start and Trio != 'UAA' and Trio != 'UAG' and Trio != 'UGA':
					Value = self.bases.index(self.sequence[i]) * 16 \
					+ self.bases.index(self.sequence[i+1]) * 4 \
					+ self.bases.index(self.sequence[i+2]) * 1
					Temp_proteins.append(self.aminoacids[Value])
					continue
			if Trio == 'UAA' or Trio == 'UAG' or Trio == 'UGA':
				self.proteins.append(Temp_proteins)
				del Temp_proteins
				Start = False
				End = True
				continue
		else:
			if Start and not End:
				#Start_index += 3 # If this is commented, the AUG will not take into
				# account for trash.
				Trash = ''
				while Start_index != i+3:
					Temp_trash = self.sequence[Start_index] + self.sequence[Start_index+1] \
					+ self.sequence[Start_index+2]
					Trash += Temp_trash
					Start_index += 3
				self.DNATrash.append(get_complement(Trash))
				self.tRNATrash.append(Trash)

		heapsort(self.frequencies)
		quicksort(self.proteins)
		#print(self.frequencies)
		#print(self.proteins)
		#print(self.tRNATrash)
		#print(self.DNATrash)

# Main program			
if __name__ == '__main__':
	DNA = DNADoble("ATGGGCAATCGGTTTGC") 
	DNA.zip()
	DNA.unzip()
	DNA.buscar("ATA")
	DNA.mitosis()
	DNA.imprimir()

	var = get_complement("ATGTTTTTCTTATTGTCTTCCTCATCGTATTACTAAATGACGATAGTAGATTGAATGTTCTAAATGTTTATGTCTTAAATGCTTAACTGAATGTTCAGCTAGATGGAGTAT")
	DNAs = DNASimple(var)
	DNAs.complement() # If you don't call this method first, the transliterate method
	transliterated = DNAs.transliterate() # will do it automatically.

	protein = tRNA(transliterated)
	protein.Translate()