#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#  
#  Copyright 2015 Rafael Caoadas <rafaelc@selje>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

from array import array

def generateReadsData(filepathReads, filepathBarcodes):

	# Get the files from disk
	print("Opening files...")
	readsFD    = open(filepathReads, 'r')
	barcodesFD = open(filepathBarcodes, 'r')

	# Count how many lines are inside
	print("Counting size file read")
	numLinesReads = sum(1 for line in readsFD)
	print(numLinesReads)
	print("Counting size file barcodes")
	numLinesBarcodes = sum(1 for line in barcodesFD)
	print(numLinesBarcodes)

	# Get the file descriptors back to the beginning
	readsFD.seek(0)
	barcodesFD.seek(0)

	
	print("Allocating memory")
	# Create the containers where we are goint to put the info for the reads
	# ID Line | Sequence | Quality
	#readsIndex     = range(numLinesReads/4)
	readsIDs       = [None] * (numLinesReads/4)
	readsSequences = [None] * (numLinesReads/4)
	readsQuality   = [None] * (numLinesReads/4)

	# We are going to use this dictionary to keep track of the IDs and which
	# barcode they use (because it is larger than the barcode string, thus
	# saving memory, we could do it the other way around too)
	readsDictionary = dict()

	# Also, we are going to keep track of every unique barcode with this set
	uniqueBarcodes = set()

	print("Filling reads data")
	# Fill the sequences data, the fastq file is divided in group of four lines
	# -- The first line is the sequence ID
	# -- The second line is the sequence
	# -- The third line is the "+" line
	# -- The forth line is the quality line
	readIDIndex = 0
	readSequenceIndex = 0
	readQualityIndex = 0
	for i in range(numLinesReads):

		nextLine = readsFD.readline()[:-1]

		# ID line
		if(i%4 == 0):
			readsIDs[readIDIndex] = nextLine
			readIDIndex += 1
			readsDictionary[nextLine.split(" ")[0]] = []

			if "21988:27580" in nextLine:
				print(nextLine)
				print(readsDictionary["@SL-MAE:A5YKC141103:A5YKC:1:1101:21988:27580"])

		# Sequence line
		if(i%4 == 1):
			readsSequences[readSequenceIndex] = nextLine
			readSequenceIndex += 1

		# Quality line
		if(i%4 == 3):
			readsQuality[readQualityIndex] = nextLine
			readQualityIndex += 1

	# Put everything into a matrix to return it later
	readsMatrix = [readsIDs, readsSequences, readsQuality]
	

	# Fill the barcode data.
	print("Filling barcodes data")
	readID = "No ID"
	for i in range(numLinesBarcodes):

		nextLine = barcodesFD.readline()[:-1] #Take out the \n

		# ID line
		if(i%4 == 0):
			readID = nextLine.split(" ")[0]

		# Sequence line
		if(i%4 == 1):
			readsDictionary[readID].append(nextLine)
			uniqueBarcodes.update(nextLine)

	return(readsDictionary,uniqueBarcodes,readMatrix)
	

def main():
	
	# TOY DATA
	forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.1.fastq"
	forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/toyBarcodes1.fastq"

	reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.2.fastq"
	reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/toyBarcodes1.fastq"

	# REAL DATA
	forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.1.fastq"
	forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.barcode_1.fastq"

	reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.2.fastq"
	reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/SN0040776/1_A5YKC.1.barcode_2.fastq"


	generateReadsData(forwardSequencesFile, forwardBarcodesFile)

	# Generate the output files

	print ("This line will be printed.")
	return 0

if __name__ == '__main__':
	main()
	

