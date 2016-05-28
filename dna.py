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
		elif nucleobase == "T" or nucleobase == 'U':
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
		print("\nBeginning search subsequence process..")
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

	var = get_compliment('ATGTTTTTCTTATTGTCTTCCTCATCGTATTACTAAATGACGATAGTAGATTGAATGTTCTAAATGTTTATGTCTTAA')
	ADNs = ADNSimple(var)
	ADNs.compliment() # If you don't call this method first, the transliterate method
	transliterated = ADNs.transliterate() # will do it automatically.

	print(transliterated)
	protein = ARNt(transliterated)
	protein.translate()