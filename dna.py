# Universidad Simón Bolívar
# CI2692: Laboratory of Algorithms and Structures II

# Project 1: Implementing of algorithms to emulate the synthesis process of 
# 			 proteins.
# Authors: David Cabeza 13-10191, Fabiola Martinez 13-10838
# Last edit: Sun 22, May 2016
import sys

def get_compliment(sequence):
	compliment = ""
	for nucleobase in sequence:
		if nucleobase == "A":
			compliment += "T"
		elif nucleobase == "T":
			compliment += "A"
		elif nucleobase == "G":
			compliment += "C"
		elif nucleobase == "C":
			compliment += "G"

	return compliment

class ADNDoble(object):
	"""Class whose instances are sequences of simple ADN. This class manipulates 
	its associated double sequence"""
	def __init__(self, simpleseq):
		self.simpleseq = simpleseq.upper()
		self.sseqlenght = len(self.simpleseq)
		self.compliment = ""
		self.doubleseq = ""

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

		self.compliment = get_compliment(self.simpleseq)
		print('The compliment for the DNA sequence given is:', self.compliment)
		print('The double sequence for the given is:', self.doubleseq)
		print('Zip process was made successfully!!')
	
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

	def mitosis(self):
		print('\nBeginning the mitosis process for the sequence', ADN.doubleseq)	
		print('Beginning to unzip the double sequence')
		seqOne = ""
		for i in range(0, len(self.doubleseq), 2):
			seqOne += self.doubleseq[i]
		seqTwo = ""
		for j in range(1, len(self.doubleseq), 2):
			seqTwo += self.doubleseq[j]
		print('Beginning to calculate the compliments')
		CseqOne = get_compliment(seqOne)
		CseqTwo = get_compliment(seqTwo)

	def buscar(self, subcadena):
		print("\nBeginning finding subsequence process..")
		compliment = get_compliment(subcadena)
		if subcadena in self.doubleseq:
			print('MATCH, we found', subcadena, 'in', self.doubleseq)
		if compliment in self.doubleseq:
			print('MATCH, we found the compliment', compliment, 'of', subcadena,'in', self.doubleseq)

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

class ADNSimple(object):
	"""Instances: Simple ADN sequence"""
	def __init__(self, sequence):
		self.sequence = sequence
		self.complim = ""

	def compliment(self):
		for nucleobase in self.sequence:
			if nucleobase == "A":
				self.complim += "T"
			elif nucleobase == "T":
				self.complim += "A"
			elif nucleobase == "G":
				self.complim += "C"
			elif nucleobase == "C":
				self.complim += "G"

	def transliterate(self):
		if not self.complim:
			self.complim == self.compliment()

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

	def translate(self):
		if 'AUG' in self.sequence:
			print('Hay un codon de inicio. Ver dónde.')

# Main program			
if __name__ == '__main__':
	ADN = ADNDoble("ATGGGCAATCGGTTTGC")
	ADN.zip()
	ADN.unzip()
	ADN.buscar("ATA")
	ADN.mitosis()
	ADN.imprimir()

	ADNs = ADNSimple("ATGGGCAATCGGTTTGC")
	ADNs.compliment() # If you don't call this method first, the transliterate method
	transliterated = ADNs.transliterate() # will do it automatically.

	print(transliterated)
	protein = ARNt(transliterated)
	protein.translate()