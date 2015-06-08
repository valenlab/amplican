import time
import operator

def getTargets(targetFilePath, firstIsForward):

	targetFD    = open(targetFilePath, 'r')

	forwardTargets = set()
	reverseTargets = set()
	comboTargets = []

	# Move the FD two lines and skip the headers
	targetFD.readline()
	targetFD.readline()

	# It could be that the first column refers to the forward barcodes or viceversa
	if(firstIsForward):
		for line in targetFD:
			columns = line.split()
			comboTargets.append(str(columns[0])+"_"+str(columns[1]))
			forwardTargets.add(columns[0])
			reverseTargets.add(columns[1])

	else:
		for line in targetFD:
			columns = line.split()
			comboTargets.append(str(columns[1])+"_"+str(columns[0]))
			forwardTargets.add(columns[1])
			reverseTargets.add(columns[0])

	return [comboTargets, forwardTargets, reverseTargets]


def main():

	start = time.time()

	#Starting string for results files
	#startString = "1_A5YKC.1."
	#startString = "1_AB7EB.1."
	#startString = "1_ABNHU.1."
	#startString = "1_AC7H3.1."
	#startString = "1_ACBKW.1."
	#startString = "1_ACCJB.1."
	startString = "1_ACBLW.1."
	

	# Target DATA
	# libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/1_library_params.txt"
	#libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0040778/1_library_params.txt"
	#libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0041333/1_library_params.txt"
	#libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0043318/1_library_params.txt"
	#libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048151/1_library_params.txt"
	#libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048152/1_library_params.txt"
	libraryFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048153/1_library_params.txt"
	
	
	# TOY DATA
	# forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/toySequences1.fastq"
	# forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/toyBarcodes1.fastq"

	# reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/toySequences2.fastq"
	# reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/toyBarcodes2.fastq"

	# REAL DATA
	# SN0040776
	# forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/1_A5YKC.1.1.fastq"
	# forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/1_A5YKC.1.barcode_1.fastq"

	# reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/1_A5YKC.1.2.fastq"
	# reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Real/SN0040776/1_A5YKC.1.barcode_2.fastq"

	# SN0040778
	# forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0040778/1_AB7EB.1.1.fastq"
	# forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0040778/1_AB7EB.1.barcode_1.fastq"

	# reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0040778/1_AB7EB.1.2.fastq"
	# reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0040778/1_AB7EB.1.barcode_2.fastq"

	# SN0041333
	#forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0041333/1_ABNHU.1.1.fastq"
	#forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0041333/1_ABNHU.1.barcode_1.fastq"

	#reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0041333/1_ABNHU.1.2.fastq"
	#reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0041333/1_ABNHU.1.barcode_2.fastq"

	# SN0043318
	#forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0043318/1_AC7H3.1.1.fastq"
	#forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0043318/1_AC7H3.1.barcode_1.fastq"

	#reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0043318/1_AC7H3.1.2.fastq"
	#reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0043318/1_AC7H3.1.barcode_2.fastq"

	# SN0048151
	#forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048151/1_ACBKW.1.1.fastq"
	#forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048151/1_ACBKW.1.barcode_1.fastq"

	#reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048151/1_ACBKW.1.2.fastq"
	#reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048151/1_ACBKW.1.barcode_2.fastq"

	# SN0048152
	#forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048152/1_ACCJB.1.1.fastq"
	#forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048152/1_ACCJB.1.barcode_1.fastq"

	#reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048152/1_ACCJB.1.2.fastq"
	#reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048152/1_ACCJB.1.barcode_2.fastq"

	# SN0048153
	forwardSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048153/1_ACBLW.1.1.fastq"
	forwardBarcodesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048153/1_ACBLW.1.barcode_1.fastq"

	reverseSequencesFile = "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048153/1_ACBLW.1.2.fastq"
	reverseBarcodesFile =  "/Home/ii/rafaelc/Desktop/CRISPRS/data/Werner/Multiples_January/SN0048153/1_ACBLW.1.barcode_2.fastq"
	

	# Prepare the file descriptor for the stats
	statsFD = open("stats.txt",'w')

	# Get the target barcodes from the library
	print("Getting target barcodes")
	targetResults = getTargets(libraryFile, False)

	comboTargets = targetResults[0]
	forwardTargets = targetResults[1]
	reverseTargets = targetResults[2]

	print("Found this many target combos:" + str(len(comboTargets)))
	print("Found this many target forwards:" + str(len(forwardTargets)))
	print("Found this many target reverses:" + str(len(reverseTargets)))

	statsFD.write("Found this many target combos:" + str(len(comboTargets)) + str("\n"))
	statsFD.write("Found this many target forwards:" + str(len(forwardTargets)) + str("\n"))
	statsFD.write("Found this many target reverses:" + str(len(reverseTargets)) + str("\n\n"))

	# Prepare the file descriptor for the combos
	barcodesDictionary = dict()
	for i in range(len(comboTargets)):
		myNewFWDFD = open(startString+comboTargets[i]+str("_unmapped.1.fastq"), 'w')
		myNewRVSFD = open(startString+comboTargets[i]+str("_unmapped.2.fastq"), 'w')

		forwardResultFD = myNewFWDFD
		reverseResultFD = myNewRVSFD

		barcodesDictionary[comboTargets[i]] = [myNewFWDFD, myNewRVSFD]

	# Get the files from disk
	print("Opening files...")
	readsFWDFD    = open(forwardSequencesFile, 'r')
	barcodesFWDFD = open(forwardBarcodesFile, 'r')
	readsRVSFD    = open(reverseSequencesFile, 'r')
	barcodesRVSFD = open(reverseBarcodesFile, 'r')

	# Count how many lines are inside
	print("Counting lines in files")
	numLinesReadsFWD = sum(1 for line in readsFWDFD)
	#numLinesReadsRVS = sum(1 for line in readsRVSFD)
	#numLinesBarcodesFWD = sum(1 for line in barcodesFWDFD)
	#numLinesBarcodesRVS = sum(1 for line in barcodesRVSFD)
	
	print("Forward Reads")
	print(numLinesReadsFWD)
	#print("Reverse Reads")
	#print(numLinesReadsRVS)
	#print("Forward Barcodes")
	#print(numLinesBarcodesFWD)
	#print("Reverse Barcodes")
	#print(numLinesBarcodesRVS)

	statsFD.write("Found this many lines:" + str(numLinesReadsFWD) + str("\n\n"))

	# Get the file descriptors back to the beginning
	readsFWDFD.seek(0)
	readsRVSFD.seek(0)
	barcodesFWDFD.seek(0)
	barcodesRVSFD.seek(0)

	# Show up all the barcodes and counting statistics
	barcodesFWD = dict()
	barcodesRVS = dict()
	for i in range(numLinesReadsFWD):
		barcodeFWDLine = barcodesFWDFD.readline()[:-1]
		barcodeRVSLine = barcodesRVSFD.readline()[:-1]

		# Sequence line; check out the barcodes
		if(i%4 == 1):

			# Increment +1 the barcodes in each or add a new one
			if(barcodeFWDLine in barcodesFWD):
				barcodesFWD[barcodeFWDLine] += 1
			else:
				barcodesFWD[barcodeFWDLine] = 1

			if(barcodeRVSLine in barcodesRVS):
				barcodesRVS[barcodeRVSLine] += 1
			else:
				barcodesRVS[barcodeRVSLine] = 1

	print("Found this many unique barcodes in Forward File")
	print(len(barcodesFWD))
	print("Found this many unique barcodes in Reverse File")
	print(len(barcodesRVS))

	statsFD.write("Found this many unique barcodes in Forward File:" + str(len(barcodesFWD)) + str("\n"))
	statsFD.write("Found this many unique barcodes in Reverse File:" + str(len(barcodesRVS)) + str("\n\n"))


	# Sort the dictionaries by value and show results
	sortedFWD = sorted(barcodesFWD.items(), key=operator.itemgetter(1))
	sortedRVS = sorted(barcodesRVS.items(), key=operator.itemgetter(1))

	showingFWD = len(forwardTargets) + 10
	showingRVS = len(reverseTargets) + 10

	print("Forward stats, showing top " +str(showingFWD) +" most representative barcodes")
	statsFD.write("Forward stats, showing top " +str(showingFWD) +" most representative barcodes" + str("\n"))
	
	i = len(barcodesFWD)-1
	targetString = ""
	while(i > len(barcodesFWD)-showingFWD-1 and i > -1):

		if(sortedFWD[i][0] in forwardTargets):
			targetString = "<- TARGET"
		else:
			targetString = ""
		
		print(str(sortedFWD[i][0]) + " = " + str(sortedFWD[i][1]*400/float(numLinesReadsFWD)) + "% " + str(targetString)) # Times 4 because the FASTQ file has one barcode for each 4 lines, *100 to display %
		statsFD.write( str(sortedFWD[i][0]) + " = " + str(sortedFWD[i][1]*400/float(numLinesReadsFWD)) + "% " + str(targetString) + str("\n"))
		
		i = i-1

	print("") #Empty line
	statsFD.write(str("\n"))

	print("Reverse stats, showing top " +str(showingRVS) +" most representative barcodes")
	statsFD.write("Reverse stats, showing top " +str(showingRVS) +" most representative barcodes" + str("\n"))
	
	i = len(barcodesRVS)-1
	targetString = ""
	while(i > len(barcodesRVS)-showingRVS-1 and i > -1):

		if(sortedRVS[i][0] in reverseTargets):
			targetString = "<- TARGET"
		else:
			targetString = ""
		
		print(str(sortedRVS[i][0]) + " = " + str(sortedRVS[i][1]*400/float(numLinesReadsFWD)) + "% " + str(targetString)) # Times 4 because the FASTQ file has one barcode for each 4 lines, *100 to display %
		statsFD.write(str(sortedRVS[i][0]) + " = " + str(sortedRVS[i][1]*400/float(numLinesReadsFWD)) + "% " + str(targetString) + str("\n"))
		
		i = i-1

	print("") #Empty line
	statsFD.write(str("\n"))

	print("Forward targets stats")
	statsFD.write("Forward targets stats\n")
	
	forwardKeys = list(forwardTargets)
	for key in forwardKeys:
		if(key in barcodesFWD):
			print(str(key) + " = " + str(barcodesFWD[key]*400/float(numLinesReadsFWD)) + "% ")
			statsFD.write(str(key) + " = " + str(barcodesFWD[key]*400/float(numLinesReadsFWD)) + "% " + str("\n"))
		else:
			print("CRITICAL ERROR!: "+key+" target barcode is not in the Forward File!")
			statsFD.write("CRITICAL ERROR!: "+key+" target barcode is not in the Forward File!" + str("\n"))

	print("") #Empty line
	statsFD.write(str("\n"))

	print("Reverse targets stats")
	statsFD.write("Reverse targets stats\n")
	
	reverseKeys = list(reverseTargets)
	for key in reverseKeys:
		if(key in barcodesRVS):
			print(str(key) + " = " + str(barcodesRVS[key]*400/float(numLinesReadsFWD)) + "% ")
			statsFD.write(str(key) + " = " + str(barcodesRVS[key]*400/float(numLinesReadsFWD)) + "% " + str("\n"))
		else:
			print("CRITICAL ERROR!: "+key+" target barcode is not in the Reverse File!")
			statsFD.write("CRITICAL ERROR!: "+key+" target barcode is not in the Reverse File!" + str("\n"))


	print("") #Empty line
	statsFD.write(str("\n"))

	# Get the file descriptors back to the beginning (again)
	readsFWDFD.seek(0)
	readsRVSFD.seek(0)
	barcodesFWDFD.seek(0)
	barcodesRVSFD.seek(0)

	# Now lets process all the files, line by line, at once.
	# Prepare some sort of crude feedback for the user
	progress = 0
	tenPercent = numLinesReadsFWD/10
	
	# In here we are going to save some information temporarely
	readFWDID = "Nothing here"
	readRVSID = "Nothing here"
	readFWDSequence = "Nothing here"
	readRVSSequence = "Nothing here"
	readFWDQuality = "Nothing here"
	readRVSQuality = "Nothing here"

	# The file descriptors where we flush the last 4 lines
	forwardResultFD = None
	reverseResultFD = None
	
	# In here we are going to save the combination of barcodes and File Descriptor for each.
	#barcodesDictionary = dict()

	skipLine = False
	currentCombo = ""
	#for i in range(0):
	for i in range(numLinesReadsFWD):

		# Update the progress
		if(i%tenPercent == 0):
			print("We are at..."+str(progress)+"%")
			progress += 10

		readFWDLine = readsFWDFD.readline()[:-1] #Take out the \n
		readRVSLine = readsRVSFD.readline()[:-1] 
		barcodeFWDLine = barcodesFWDFD.readline()[:-1]
		barcodeRVSLine = barcodesRVSFD.readline()[:-1]

		# Figure it out in which line are we:
		# ID line; we need to save the content for later
		if(i%4 == 0):
			readFWDID = readFWDLine
			readRVSID = readRVSLine
			
		# Sequence line; check out the barcodes and save the content
		if(i%4 == 1):

			# Save the sequence
			readFWDSequence = readFWDLine
			readRVSSequence = readRVSLine

			# Get the barcodes and figure it out the ID (combo)
			forwardBarcode = barcodeFWDLine
			reverseBarcode = barcodeRVSLine
			currentCombo = forwardBarcode + str("_") + reverseBarcode

			# Do we have a file descriptor for this file already?
			if(currentCombo in barcodesDictionary):
				# If so, retrieve it, and prepare to write into it
				forwardResultFD = barcodesDictionary[currentCombo][0]
				reverseResultFD = barcodesDictionary[currentCombo][1]

				skipLine = False

			else:
				skipLine = True
				
			
		# Quality line; we flush everything
		if(i%4 == 3):

			if(skipLine == False):

				# Save the quality
				readFWDQuality = readFWDLine
				readRVSQuality = readRVSLine

				# Write the stuff into the files
				forwardResultFD.write(readFWDID + str("\n"))
				forwardResultFD.write(readFWDSequence + str("\n"))
				forwardResultFD.write("+" + str("\n"))
				forwardResultFD.write(readFWDQuality + str("\n"))

				reverseResultFD.write(readRVSID + str("\n"))
				reverseResultFD.write(readRVSSequence + str("\n"))
				reverseResultFD.write("+" + str("\n"))
				reverseResultFD.write(readRVSQuality + str("\n")) 

	# Now close all the opened FDs
	allFDKeys = barcodesDictionary.keys()
	for key in allFDKeys:
		barcodesDictionary[key][0].close()
		barcodesDictionary[key][1].close()
	
	statsFD.close()


	end = time.time()
	print("Total time")
	print (end - start)
	
	return 0

if __name__ == '__main__':
	main()
	
