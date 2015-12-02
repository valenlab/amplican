# This function takes a plot folder and generate the proper HTML
# documents. It can also copy the documents to a given folder inside
# your server.
#
# Parameters:
#
# (string) PLOTFOLDER     The path to your plot folder. For example:
#
#                         /Home/johndoe/.../results/run_XYZ/Analysis_XYZ/Plots_XYZ
#
#                         Note that the plot folder is inside the
#                         analysis folder, which is inside the
#                         alignments folder. So if you have different
#                         plots, and different analysis for a given
#                         sample you can generate easily an HTML for
#                         each.
#
# (string) SERVERFOLDER   The path to your plot server folder, where the
#                         documents are going to get copied.
#                         For example:
#
#                         "/net/public_html/johndoe/myRun003"
#
#                         If SERVERFOLDER is equal to "" (default), it
#                         won't try to copy anything and leave the
#                         documents inside the HTML folder, which will
#                         be inside PLOTFOLDER.
#
# (string) SERVERURL      Where are going to be located your results
#                         if you want to access them via web browser.
#                         For example:
#
#                         "http://www.myserver.no/~johndoe/myRun003"
#
# (bool) STAMP            We can generate an HTML with a time stamp if
#                         so we desire like with the rest of the
#                         program. Generally, only one is a good idea.
#                         If STAMP == FALSE only one folder will be
#                         created.
#

ampliCanHTML <- function(PLOTFOLDER, SERVERURL = "", SERVERFOLDER = "",
                         STAMP = FALSE){


	############################################
	#Summons the libraries and get to use functions from there.
	############################################
	{
	  #library("png") # HTML functions will use this to resize RAW pngs into tiny previews
	  library("raster")
	  source("libraries/configHTML.R") # Function that deals with the config files and get running everything.
	}

	############################################
	# Init some constants
	############################################
	{

  #SERVERURL       <- "http://www.ii.uib.no/~eivindv/Jason2"
  #SERVERURL       <- "http://www.ii.uib.no/~eivindv/JasonAfterChristmas"
  #SERVERURL       <- "http://www.ii.uib.no/~eivindv/JasonMay"
  

  #SERVERURL       <- "http://www.ii.uib.no/~eivindv/Danielle_F"
  #SERVERURL       <- "http://www.ii.uib.no/~eivindv/Danielle_R"
  
	  # Get also the CSS file with all the cool graphical stuff
	  CSSPATH <- "/../../res/HTML/style.css"
		
	  # Get also the JS file with all the cool functionality
	  JSPATH  <- "/../../res/HTML/scripts.js"
	  
	  # Several resources for the document
	  LOGO    <- "/../../res/HTML/Logo.png"

	  # Find out the folders
	  
	  # -- DATA FOR THE PLOTS
	  DATAFOLDER <- paste(PLOTFOLDER,"/Data",sep='',collapse='')
	  
	  # -- ANALYSIS
	  plotTree            <- strsplit(PLOTFOLDER, "/", fixed = TRUE)[[1]]
	  plotTreeDeep        <- length(plotTree)
	  analysisTree        <- plotTree[1:plotTreeDeep -1]
	  ANALYSISFOLDER      <- paste(analysisTree, collapse='/')
	  
	  # -- ALIGNMENTFOLDER and RESULTSFOLDER
	  analysisTreeDeep    <- length(analysisTree)
	  alignmentTree       <- analysisTree[1:analysisTreeDeep -1]
	  RESULTSFOLDER       <- paste(alignmentTree, collapse='/')
	  ALIGNMENTFOLDER     <- list.files(path = RESULTSFOLDER, pattern = "_alignments", full.names = TRUE)
	  ALIGNMENTFOLDERNAME <- list.files(path = RESULTSFOLDER, pattern = "_alignments", full.names = FALSE)
	  
	  #print(ANALYSISFOLDER)
	  #print(ALIGNMENTFOLDER)
	  #print(RESULTSFOLDER)
	  
	  # TODO:
	  # Get the format
	  FORMAT <- ".png"
	  # Get the group by info
	  GROUPBY <- 20
	  LINKBY  <- 5
	}

	############################################
	# Algorith starts here
	############################################
	{

	  #------------------------------------------------------
	  # CREATE THE FOLDER STRUCTURE
	  #------------------------------------------------------
	  {
		# Now, we "construct" the server folder.
		# It will have this format:
		# /Plot Folder
		# ---- /Plots withe the cutrates
		# ---- /Plots with the barcodes
		# ---- /Plots with frameshifts
		# ---- /Plots with the archplots and histograms
		# /Data Folder with all the TXT with the RAW data of the plots
		# Big tar.gz with all the alignments
		# Big tar.gz with all the analysis
		#/HTML Folder
		# ---- index.html
		# ---- barcodes.html
		# ---- Results1.html
		# ---- ...
		# ---- ResultsX.html
		
		#Create the Server main folder
		dir.create(file.path(SERVERFOLDER), showWarnings = FALSE)
		
		# Create a folder called Plot Folder in the Server Folder
		plotsPath <- paste(SERVERFOLDER, "/Plots", sep='', collapse='')
		dir.create(file.path(plotsPath), showWarnings = FALSE)
		
		# Inside that folder create the 4 subfolders for the plots
		plotsCutratesPath   <- paste(plotsPath, "/Cutrates",   sep='',collapse='')
		plotsBarcodesPath   <- paste(plotsPath, "/Barcodes",   sep='',collapse='')
		plotsFrameshiftPath <- paste(plotsPath, "/Frameshift", sep='',collapse='')
		plotsArchplotsPath  <- paste(plotsPath, "/Archplots",  sep='',collapse='')
		
		plotsCutratesURL   <- paste(SERVERURL, "/Plots/Cutrates",   sep='',collapse='')
		plotsBarcodesURL   <- paste(SERVERURL, "/Plots/Barcodes",   sep='',collapse='')
		plotsFrameshiftURL <- paste(SERVERURL, "/Plots/Frameshift", sep='',collapse='')
		plotsArchplotsURL  <- paste(SERVERURL, "/Plots/Archplots",  sep='',collapse='')
		
		dir.create(file.path(plotsCutratesPath), showWarnings = FALSE)
		dir.create(file.path(plotsBarcodesPath), showWarnings = FALSE)
		dir.create(file.path(plotsFrameshiftPath), showWarnings = FALSE)
		dir.create(file.path(plotsArchplotsPath), showWarnings = FALSE)
		
		# Create a folder called Data in the Server Folder
		dataPath  <- paste(SERVERFOLDER, "/Data", sep='', collapse='')
		dataURL  <- paste(SERVERURL, "/Data", sep='', collapse='')
		
		dir.create(file.path(dataPath), showWarnings = FALSE)
	  
	  }

	  #------------------------------------------------------
	  # COPY THE FILES INTO THE PLOT FOLDERS
	  #------------------------------------------------------
	  {
		allFiles <- list.files(path = PLOTFOLDER, pattern = FORMAT, full.names = TRUE)
		
		# For all files, see where they go
		for(i in 1:length(allFiles)){
		  
		  # Copy the cut rates  
		  if(length(grep("Cut_Rates_", allFiles[i]) == 1)){
			file.copy(from=allFiles[i], to=plotsCutratesPath, overwrite = TRUE, copy.mode = TRUE)
		  }
		  else{
		  
			# Copy the frameshift
			if(length(grep("Frameshift_", allFiles[i]) == 1)){
			  file.copy(from=allFiles[i], to=plotsFrameshiftPath, overwrite = TRUE, copy.mode = TRUE)
			}
			
			else{
			  # Copy the barcodes plots, which is either barcodes for assigned, filters, or unnasigned
			  if(length(grep("barcode_", allFiles[i]) == 1) || length(grep("_unnasigned.", allFiles[i]) == 1)){
				file.copy(from=allFiles[i], to=plotsBarcodesPath, overwrite = TRUE, copy.mode = TRUE)
			  }
			  else{
				
				# Copy the Gene ID specific plots (archplots, histograms, etc...)
				file.copy(from=allFiles[i], to=plotsArchplotsPath, overwrite = TRUE, copy.mode = TRUE)
			  
			  }
			
			}
		  
		  }
		  
		}
		
		dataFiles <- list.files(path = DATAFOLDER, full.names = TRUE)
		for(i in 1:length(dataFiles)){
		  file.copy(from=dataFiles[i], to=dataPath, overwrite = TRUE, copy.mode = TRUE)
		}
		
	  
	  }
	  
	  #------------------------------------------------------
	  # CREATE THE TAR.GZ FOR THE ALIGNMENTS AND THE ANALYSIS
	  # -- This will also copy the files into the server folder
	  #------------------------------------------------------
	  {
	   
		#TODO: Fix this
		
		# We need to change the working folder to the folder with the alignments.
		# Otherwise we get a zip with the root set as the path of the current working directory
		# So the zip files contains a folder call /home/ii/.../alignments instead of just /alignments for example
		current_folder <- getwd()
		
		setwd(RESULTSFOLDER)
		allFiles <- list.files(path = RESULTSFOLDER, pattern = "", full.names = FALSE)
		tarName  <- paste(plotTree[plotTreeDeep-2],".tar.gz",collapse='',sep='')
		tar(tarName, files = allFiles, compression = "gzip", extra_flags = "--format=ustar")
		tarFilePath <- paste(RESULTSFOLDER,"/",tarName,collapse='',sep='')
		file.copy(from=tarFilePath, to=SERVERFOLDER, overwrite = TRUE, copy.mode = TRUE)
		
		# Set the working directory to the alignments, make the gz and move it to the server folder
	#     setwd(ALIGNMENTFOLDER)
	#     tar("alignments.tar.gz", files = ALIGNMENTFOLDERNAME, compression = "gzip")
	#     alignmentFile <- paste(RESULTSFOLDER,"/alignments.tar.gz", sep='', collapse='')
	#     file.copy(from=alignmentFile, to=SERVERFOLDER, overwrite = TRUE, copy.mode = TRUE)
	#     
	#     # Set the working directory to the analysis, make the gz and move it to the server folder
	#     setwd(ANALYSISFOLDER)
	#     tar("analysis.tar.gz", files = ALIGNMENTFOLDERNAME, compression = "gzip")
	#     analysisFile <- paste(RESULTSFOLDER,"/analysis.tar.gz", sep='', collapse='')
	#     file.copy(from=analysisFile, to=SERVERFOLDER, overwrite = TRUE, copy.mode = TRUE)
		
		
		# Set the working directory to the analysis
		# THIS REGEX DO NOT WORK!
		# list.files(path = ANALYSISFOLDER, pattern = "^((?!Plots_).)*$", full.names = TRUE)
		#plotFolders <- list.files(path = ANALYSISFOLDER, pattern = "Plots_", full.names = TRUE)
		
		
		
		# Revert changes in the working directory
		setwd(current_folder)
		
	  }
	  
	  #------------------------------------------------------
	  # COPY THE BARCODE TABLE INTO THE DATA FOLDER
	  # TODO: Copy this with everything else?
	  #------------------------------------------------------
	  {
		
		barcodeFile <- paste(ALIGNMENTFOLDER, "/barcodeFile_results.txt", sep='', collapse='')
		file.copy(from=barcodeFile, to=dataPath, overwrite = TRUE, copy.mode = TRUE)
		
	  }
	  
	  #------------------------------------------------------
	  # GENERATE THE HTML CODE FOR EACH PAGE
	  #------------------------------------------------------
	  {

		# Create the folder where the .htmls is going to be
		if(STAMP == TRUE){
		  currentTime <- Sys.time()
		  timeStamp   <-  strftime(currentTime,"%Y%m%d%H%M%S")
		  htmlPath    <- paste(PLOTFOLDER, "/HTML_",timeStamp, sep = '')
		  dir.create(file.path(htmlPath), showWarnings = TRUE)
		}else{
		  htmlPath    <- paste(PLOTFOLDER, "/HTML", sep = '')
		  dir.create(file.path(htmlPath), showWarnings = FALSE)
		}
		
		# Copy the CSS and JS files in there
		file.copy(from=CSSPATH, to=htmlPath, copy.mode = TRUE, overwrite = TRUE)
		file.copy(from=JSPATH,  to=htmlPath, copy.mode = TRUE, overwrite = TRUE)
		
		# Copy any other resource needed
		file.copy(from=LOGO, to=htmlPath, copy.mode = TRUE, overwrite = TRUE)
		
		# Figure it out how many groups we have
		# -- Get the config file from the analysis folder
		# -- Count how many rows we have
		analysisDF <- read.table(paste(ANALYSISFOLDER,  "/analysis_results.txt",     sep=''), header=TRUE, sep="\t")
		barcodeDF  <- read.table(paste(ALIGNMENTFOLDER, "/barcodeFile_results.txt" , sep=''), header=TRUE, sep="\t")
		
		# We don't need the timing stuff in here, just keep the interesting columns
		keeps <- c("ID","Barcode","Forward_Primer","Reverse_Primer","Genome","Sum_Is_Sequence","Cut_Rate","Forward_No_Deletions","Forward_In_Frameshift","Forward_Out_Frameshift","Forward_Multiple_Deletions","Reverse_No_Deletions","Reverse_In_Frameshift","Reverse_Out_Frameshift","Reverse_Multiple_Deletions","Sum_Forward_Cuts","Sum_Reverse_Cuts","Sum_Is_Read","Sum_Unique_Cuts_Passed_ALLIN","Sum_Unique_Cuts_Passed_SAMECUT","Sum_Unique_Cuts_Passed_PRIMERDIMER","Sum_Unique_Reads_Passed_ALLIN","Sum_Unique_Reads_Passed_SAMECUT","Sum_Unique_Reads_Passed_PRIMERDIMER","Sum_Unique_Cuts_After_ALLIN","Sum_Unique_Cuts_After_SAMECUT","Sum_Unique_Cuts_After_PRIMERDIMER","Sum_Unique_Reads_After_ALLIN","Sum_Unique_Reads_After_SAMECUT","Sum_Unique_Reads_After_PRIMERDIMER")
		anaylisLittleDF <- analysisDF[keeps]
		
		# Sort the dataframe at the same order we got the original config file
		# This could be that the GeneID is not sorted alphabetically, we don't care about it.
		anaylisLittleDF <- anaylisLittleDF[ order(as.numeric(row.names(anaylisLittleDF))), ]
		barcodeDF       <- barcodeDF[ order(as.numeric(row.names(barcodeDF))), ]
		
		# Divide the stuff
		analyisisSubDF <- split(anaylisLittleDF, factor(sort(rank(row.names(anaylisLittleDF))%%ceiling(nrow(anaylisLittleDF)/GROUPBY))))
		barcodeSubDF   <- split(barcodeDF,       factor(sort(rank(row.names(barcodeDF))      %%ceiling(nrow(barcodeDF)      /GROUPBY))))
		
		# -----------------
		# INDEX
		# -----------------
		print("Making HTML Index")
		makeIndexHTML3(htmlPath, analyisisSubDF, SERVERURL)
		
		# -----------------
		# STATISTICS
		# -----------------
		print("Making HTML Statistics")
		print(SERVERURL)
		makeStatisticsHTML2(htmlPath, PLOTFOLDER, ANALYSISFOLDER, ALIGNMENTFOLDER, RESULTSFOLDER, SERVERURL)
		
		# -----------------
		# BARCODES
		# -----------------
		print("Making HTML Barcode")
		makeBarcodeHTML2(htmlPath, ANALYSISFOLDER, PLOTFOLDER, barcodeSubDF, plotsBarcodesURL, dataURL, SERVERURL, LINKBY)
			
		# -----------------
		# Analysis
		# -----------------
		print("Analysis")    
		makeCutRatesHTML3(htmlPath, ANALYSISFOLDER, PLOTFOLDER, LINKBY, analyisisSubDF,
						  plotsCutratesURL, plotsFrameshiftURL, dataURL, SERVERURL)
	   
		# -----------------
		# Results
		# -----------------
		print("Results")    
		makeResultsHTML4(htmlPath, ANALYSISFOLDER, plotsArchplotsPath, analyisisSubDF, plotsArchplotsURL, dataURL, SERVERURL)
		
	  }
	  
	  #------------------------------------------------------
	  # MOVE THE HTML FOLDER INTO THE SERVER FOLDER
	  #------------------------------------------------------
	  {
	  
		if(SERVERFOLDER != ""){

			htmlServerPath <- paste(SERVERFOLDER, "/HTML", sep='', collapse='')
			dir.create(file.path(htmlServerPath), showWarnings = FALSE)
			
			allHTMLFiles <- list.files(path = htmlPath, full.names = TRUE)
			
			for(i in 1:length(allHTMLFiles)){
			  file.copy(from=allHTMLFiles[i], to=htmlServerPath, copy.mode = TRUE, overwrite = TRUE)
			  
			}		

		
		}
		
	  }
	  
	}

  
}







