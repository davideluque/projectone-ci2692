# Universidad Simón Bolívar
# CI2692: Laboratory of Algorithms and Structures II
# Project 1: Implementing of algorithms to emulate the synthesis process of 
# 			 proteins.
# Authors: David Cabeza 13-10191, Fabiola Martinez 13-10838
# Last edit: Mon 06, May 2016

import sys
from random import randint
from colorama import *
init()

################################################################################
#								METHODS										   #
################################################################################

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

######################### transliterate_seq METHOD #############################
# Description: This method receives a RNA sequence and transliterates it by
#			   changing Uracyl for Timine and viceversa
# Input: Receives a sequence of DNA to be complemented 
# Output: Returns the complement of the given DNA sequence
################################################################################
def transliterate_seq(sequence):
	temp = ""
	for base in sequence:
		if base == "U": temp += "T"
		elif base == "T": temp += "U"
		else: temp += base

	return temp

############################# read_sequence METHOD #################################
# Description: Reads a file that contains sequences of DNA and stores them in an
#			   array.
# Input: Receives file name.
# Output: An array with the sequences.
################################################################################
def read_sequence(DnaFile):
	
	DNA_sequences = []
	
	with open(DnaFile, 'r') as f:

		for line in f:
			line = line.rstrip()
			if line:
				DNA_sequences.append(line)
		
		f.close()

	return DNA_sequences

def loop_heap(heap):
	for heaps in heap:
		print('Codon:', heaps[0], 'with frequency:', heaps[1])

########################### quicksort METHOD ###################################
# Description: This quicksort sorts the array given by the lenght of the arrays 
#			   If arrays have the same lenght, the sorting is made by alphabetic
#			   order.
# Input: Array to sort, start index(optional), end index(optional)
# Output: Array sorted decreasingly.
################################################################################
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

def mergesort_merge(A, p, q, r):
	n = q - p + 1
	m = r - q 

	L = [A[p+i-1] for i in range(1, n+1)]
	R = [A[q+i] for i in range(1, m+1)]
	
	i = j = 0

	for k in range(p, r+1):
		if i == len(L): A[k] = R[j]
		elif j == len(R): A[k] = L[i]
		elif L[i] <= R[j]:
			A[k] = L[i]
			i += 1
		else:
			A[k] = R[j]
			j += 1

def mergesort(A, p = None, r = None):
	p = p if p else 0 
	r = r if r else len(A) - 1

	if p < r:
		q = (p + r) // 2
		mergesort(A, p, q)
		mergesort(A, q+1, r)
		mergesort_merge(A, p, q, r)

	return A

################################################################################
#								CLASSES										   #
################################################################################

class DNADouble(object):
	"""Class whose instances are simple DNA sequences. Manipulates associated 
	double sequence"""
	def __init__(self, sequence):
		self.sequence = sequence.upper()
		self.lenght = len(self.sequence)
		self.complement = ""
		self.double = ""

	############################# zip METHOD ###################################
	# Description: Receives a simple DNA sequence, calculates its complement and
	#			   returns the double DNA sequence.
	# Input: Simple DNA sequence
	# Output: Double DNA sequence
	############################################################################
	def zip(self):
		print('\nBeginning to zip the sequence', DNA.sequence)
		
		for nucleobase in self.sequence:
			
			if nucleobase == "A":
				complement = "T"
				self.complement += complement
				self.double += nucleobase + complement
			
			elif nucleobase == "T":
				complement = "A"			
				self.complement += complement
				self.double += nucleobase + complement
			
			elif nucleobase == "G":
				complement = "C"			
				self.complement += complement
				self.double += nucleobase + complement
			
			elif nucleobase == "C":
				complement = "G"		
				self.complement += complement
				self.double += nucleobase + complement

		print('Complement:', self.complement)
		print('Double sequence', self.double)

		return self.double
	
	############################# unzip METHOD #################################
	# Description: Receives a double DNA sequence and splits it into two simple
	#			   sequences
	# Input: Double DNA sequence
	# Output: Simple DNA sequences
	############################################################################
	def unzip(self):
		if not self.double:
			print('Couldn\'t unzip sequence. Please call first zip method.')
			return None

		else:
			print('\nBeginning to unzip the sequence', self.double)	
			
			simple_a = ""
			for i in range(0, len(self.double), 2):
				simple_a += self.double[i]
			
			simple_b = ""
			for j in range(1, len(self.double), 2):
				simple_b += self.double[j]
			
			print('Simple sequences are:', simple_a, 'and', simple_b)

			return simple_a, simple_b

	############################# mitosis METHOD ###############################
	# Description: Creates a copy of the double sequence through mitosis
	# Input: Double DNA sequence
	# Output: Simple DNA sequences
	############################################################################
	def mitosis(self):
		print('\nBeginning mitosis for:', self.double)
		simple_a, simple_b = self.unzip()	

		simple_a_complement = get_complement(simple_a)
		simple_b_complement = get_complement(simple_b)

		return simple_a_complement, simple_b_complement
	
	############################# search METHOD ###############################
	# Description: Receives a simple DNA sequence and looks for a match in DNA.
	#			   The search is made simultaneously for the sequence and its
	#			   complement
	# Input: Simple sequence
	# Output: Returns True if a match exists
	############################################################################
	def search(self, subsequence):
		print("\nSearching", subsequence, 'in', self.double, '...')
		
		i = 0
		for j in range(0, len(self.double)):
			if i != len(subsequence):
				if subsequence[i] == self.double[j]:
					Match = True
					i += 1
					continue
				else:
					Match = False
					i = 0
					continue
			else:
				match_index = j - len(subsequence)
				print('Match, the subsequence begins in position', match_index, \
				'of the double sequence.')
				return True 
				break
		else:
			if not Match:
				print('Subsequence was not found')
				return False

		complement = get_complement(subsequence)
		if complement in self.double:
			print('MATCH, we found the complement', complement, 'of', subsequence)

	############################# printpairs METHOD ############################
	# Description: Prints double DNA sequence in pairs
	# Input: Double DNA sequence
	# Output: Pairs in parentheses () of double DNA
	############################################################################
	def printpairs(self):
		print()
		for i in range(0, len(self.double), 2):
			print('('+ Fore.RED + self.double[i], Fore.BLUE + self.double[i+1], Style.RESET_ALL + ')', sep ='', end='')
		print()

	############################### write METHOD ###############################
	# Description: Add (creates) to a (a) file with the double DNA sequence
	# Input: Name of the file
	# Output: True if the sequence is written 
	############################################################################
	def write(self, file):
		with open(file, 'a') as f:
			f.write(self.double + '\n')
			f.close()		
		assert(f.closed)

		return True

class DNASimple(object):
	"""Instances: Simple DNA sequence"""
	def __init__(self, sequence):
		self.sequence = sequence
		self.sequence_complement = ""

	######################## complement METHOD #################################
	# Description: This method receives a DNA sequence and complements it by
	#			   changing Timine for Adenine, Guanine for Cytosine and vice 
	#			   versa.
	# Input: Receives a sequence of DNA to be complemented 
	# Output: Returns the complement of the given DNA sequence
	############################################################################
	def complement(self):
		
		for nucleobase in self.sequence:

			if nucleobase == "A":
				temp_complement = "T"
				self.sequence_complement += temp_complement
			
			elif nucleobase == "T":
				temp_complement = "A"
				self.sequence_complement += temp_complement
			
			elif nucleobase == "G":
				temp_complement = "C"
				self.sequence_complement += temp_complement
			
			elif nucleobase == "C":
				temp_complement = "G"
				self.sequence_complement += temp_complement

		return self.sequence_complement

	######################## transliterate METHOD ##############################
	# Description: This method changes Timine by Uracil in the complement of the
	#			   simple sequence.
	# Input: Receives the complement of a simple sequence
	# Output: Returns the sequence as tRNA with Timine changed by Uracil
	############################################################################
	def transliterate(self):
		if not self.sequence_complement:
			self.sequence_complement == self.complement()

		self.tRNA = ""
		for nucleobase in self.sequence_complement:
			if nucleobase == "T":
				self.tRNA += "U"
			
			else:
				self.tRNA += nucleobase

		return self.tRNA

	############################### write METHOD ###############################
	# Description: Add (creates) to a (a) file with the simple DNA sequence
	# Input: Name of the file
	# Output: True if the sequence is written
	############################################################################
	def write(self, file):
		with open(file, 'a') as f:
			f.write(self.sequence + '\n')
			f.close()		
		assert(f.closed)

		return True

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
			if Start and (Trio == 'UAA' or Trio == 'UAG' or Trio == 'UGA'):
				self.proteins.append(Temp_proteins)
				del Temp_proteins
				Start = False
				End = True
				continue
		else:
			if Start and not End:
				#Start_index += 3 # If this is commented, the AUG will be taken into
				# account for trash.
				Trash = ''
				while Start_index != i+3:
					Temp_trash = self.sequence[Start_index] + self.sequence[Start_index+1] \
					+ self.sequence[Start_index+2]
					Trash += Temp_trash
					Start_index += 3
				self.DNATrash.append(transliterate_seq(Trash))
				self.tRNATrash.append(Trash)

		heapsort(self.frequencies)
		quicksort(self.proteins)
		#print(self.DNATrash)
		mergesort(self.DNATrash)
		#print(self.DNATrash)
		#print(self.frequencies)
		#print(self.proteins)
		#print(self.tRNATrash)
		loop_heap(self.frequencies)

	############################### write METHOD ###############################
	# Description: Add (creates) to a (a) file with the array of proteins
	# Input: Name of the file
	# Output: True if the proteins are written
	############################################################################
	def write(self, file):
		with open(file, 'a') as f:
			f.write(self.proteins + '\n')
			f.close()		
		assert(f.closed)

		return True

################################################################################
#								MAIN PROGRAM								   #
################################################################################			
if __name__ == '__main__':
	DNA = DNADouble("ATGGGCAATCGGTTTGC") 
	DNA.zip()
	DNA.unzip()
	DNA.search("ATA")
	DNA.search("TGA")
	DNA.mitosis()
	DNA.printpairs()
	DNA.write('cadenas.txt')	
	
	'''
	aminoacidos = read_sequence('complejo.txt')
	for i in aminoacidos:
		DNAs = DNASimple(get_complement(i))
		DNAs.complement()
		transliterated = DNAs.transliterate()
		protein = tRNA(transliterated)
		protein.Translate()
	'''

	#print('y el de siempre...')		
	var = get_complement("ATGTTTTTCTTATTGTCTTCCTCATCGTATTACTAAATGACGATAGTAGATTGAATGTTCTAAATGTTTATGTCTTAAATGCTTAACTGAATGTTCAGCTAGATGGAGTAT")
	DNAs = DNASimple(var)
	DNAs.complement() # If you don't call this method first, the transliterate method
	transliterated = DNAs.transliterate() # will do it automatically.

	protein = tRNA(transliterated)
	protein.Translate()
