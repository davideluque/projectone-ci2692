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
	"""Clase que tiene como instancias cadenas de ADN simple y manipula su 
	cadena de ADN doble asociada"""
	def __init__(self, simpleseq):
		self.simpleseq = simpleseq.upper()
		self.sseqlenght = len(self.simpleseq)
		self.compliment = ""
		self.doubleseq = ""
		self.dseqlenght = len(self.doubleseq)

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
		print('\nBeginning the mitosis proccess for the sequence', ADN.doubleseq)	
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


ADN = ADNDoble("ATGGGCAATCGGTTTGC")

ADN.zip()
ADN.unzip()
ADN.buscar("ATA")
ADN.mitosis()
ADN.imprimir()