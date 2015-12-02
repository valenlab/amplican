# First, we need to create the index HTML.
# In here, we are going to have this:
# --- Link for the statistics results
# --- Link for the barcodes results
# --- --- Link for the assigned plot
# --- --- Link for the filters plot
# --- Link for the config table result
# --- --- Link for the cut rate plot
# --- --- Link for the frameshift plot
# --- --- Results 1 From AS392 to G4930
# --- --- Results 2 From G4931 to T2039
#         ....
# --- --- Results N From X     to Z
makeIndexHTML3 <- function (htmlPath, analyisisSubDF, SERVERURL){
  
  # Create the HTML source file. We are going to write here the HTML result for this function.
  # at the end of this function, this will contain the entire html
  sourceFD <-  file(paste(htmlPath,"/index.html",collapse = '',sep = '') , open="w")
  
  # Figure it out the name for the tar file
  htmlTree            <- strsplit(htmlPath, "/", fixed = TRUE)[[1]]
  htmlTreeDeep        <- length(htmlTree)
  tarName             <- htmlTree[htmlTreeDeep-3]
  tarName             <- paste(tarName,".tar.gz",collapse='',sep='')
  
  htmlSource <-paste(" <!DOCTYPE html>
                     <html>
                     
                     <head>
                     
                     <link rel='stylesheet' type='text/css' href='style.css'>
                     
                     </head>
                     
                     <body>
                     
                     <div id='indexBodyContainer' class='bodyContainer'>
                     
                     <div id='indexHeader' class='header'>
                     
                     <div id='indexLogo' class='logo'>
                     
                     <img src='Logo.png' alt='I am trap in the server room, please send help!'> 
                     
                     </div>
                     
                     </div>
                     
                     <div id='indexRest' class='rest'>
                     
                     <div id='index' class='map'>
                     
                     <div id='indexBackToParent' class='parentTitle boxTitle'>
                     <div id='linkIndexTitle' class='boxElement'>
                     <a href='",SERVERURL,"'> Parent Folder </a>
                     </div>
                     </div>
                     
                     <div id='indexTitle' class='indexTitle boxTitle'>
                     <div id='linkIndexTitle' class='boxElement'>
                     Index
                     </div>
                     </div>
                     
                     </div>
                     
                     <div id='indexTableOfContent' class='TOC'>
                     
                     <div id='indexTOCTitle' class='TOCTitle'>
                     <h2>TABLE OF CONTENT</h2>
                     </div>
                     
                     <div id=indexTOCContent class='TOCContent'>
                     
                     <ul>
                     <li class='indexLinkLevel_1'> <a href='statistics.html'> <span class='indexNumber'>1</span> <span class='indexText'>Statistics</span> </a> </li>
                     <li class='indexLinkLevel_1'> <a href='barcodes.html'>   <span class='indexNumber'>2</span> <span class='indexText'>Barcodes  </span> </a> 
                     <ul>
                     <li class='indexLinkLevel_2'> <a href='barcodes.html#BarcodesDefinitions'>  <span class='indexNumber'>2.1</span> <span class='indexText'>Definitions</span>    </a></li>
                     <li class='indexLinkLevel_2'> <a href='barcodes.html#BarcodesPlots'> <span class='indexNumber'>2.2</span> <span class='indexText'>Plots</span> </a></li>
                     <li class='indexLinkLevel_2'> <a href='",SERVERURL,"/barcodeFile_results.txt'>        <span class='indexNumber'>2.4</span> <span class='indexText'>Barcode info (TXT)</span> </a></li>
                     </ul> </li>
                     
                     <li class='indexLinkLevel_1'> <a href='frameshift.html'>  <span class='indexNumber'>3</span> <span class='indexText'>Analysis</span> </a>
                     <ul>
                     <li class='indexLinkLevel_2'> <a href='frameshift.html#AnalysisDefinitions'> <span class='indexNumber'>3.1</span> <span class='indexText'> Definitions </span> </a></li>
                     <li class='indexLinkLevel_2'> <a href='frameshift.html'>                     <span class='indexNumber'>3.2</span> <span class='indexText'> Plots       </span> </a></li>
                     <li class='indexLinkLevel_2'> <a href='",SERVERURL,"/",tarName,"'>       <span class='indexNumber'>3.3</span> <span class='indexText'> All data (TAR.GZ) </span> </a></li>
                     
                     </ul> </li>

                     <li class='indexLinkLevel_1'> <a href='#Archplots'> <span class='indexNumber'> 4 </span> <span class='indexText'>Archplots</span> </a>
                     <ul>
                     
                     ", sep='', collapse='')
  
  writeLines(htmlSource, sourceFD)
  # Fill the results index
  
  # We need to adjust the names with spaces so the index looks pretty, lets find maximum of characters we need
  maxLengthName <- 0
  for (i in 1:length(analyisisSubDF)){
    
    if(nchar(as.character(analyisisSubDF[[i]][1,]$"ID")) > maxLengthName){
      maxLengthName <- nchar(as.character(analyisisSubDF[[i]][1,]$"ID"))
    }
    
  }
  
  # We have many groups, for each group, put it in the analysis index
  for (i in 1:length(analyisisSubDF)){
    
    startLength <- nchar(as.character(analyisisSubDF[[i]][1,]$"ID"))
    extraString <- paste(rep("_",(maxLengthName - startLength)),collapse='',sep='')
    
    start <- paste(as.character(analyisisSubDF[[i]][1,]$"ID"),extraString,collapse='',sep='')
    end   <- as.character(analyisisSubDF[[i]][nrow(analyisisSubDF[[i]]), ]$"ID")
    
    htmlSource <-paste("<li class='indexLinkLevel_3'>
                       <a href= 'htmlResults_",i,".html'>
                       <span class='indexText'>", start ," ----> ", end ,"</span>
                       </a></li>",collapse='',sep='')
    
    writeLines(htmlSource, sourceFD)
    
  }
  htmlSource <-" </ul> </li>  "
  
  writeLines(htmlSource, sourceFD)
  
  
  
  # Close the:
  # -- TOCContent div
  # -- TOC div
  # -- rest div
  # -- body container div
  # -- body html
  # -- html
  
  htmlSource <-paste(" 
                              </div>

                            </div>

                          </div>

                        </div>
                     
                      </body>
                    </html>
                     
                     ")
  
  writeLines(htmlSource, sourceFD)
  
  
  
  # End the BODY and HTML tags
  htmlSource <- "</body> </html>"
  writeLines(htmlSource, sourceFD)
  close(sourceFD)
  
}

# Get the constant for a particular alignment run
# It will return in this order: 
# totalProcessors, skipBadNucleotides, scoringMatrix, gapOpening, gapExtension
getAlignmentConstant <- function(AlignmentFolder){

  textFile <- readLines(paste(AlignmentFolder,"/MasterLog.txt", collapse='', sep=''),encoding="UTF-8")
  
  totalProcessors    <- as.numeric(strsplit(textFile[1], " ", fixed = TRUE)[[1]][9])
  skipBadNucleotides <- strsplit(textFile[2], " ", fixed = TRUE)[[1]][6]
  scoringMatrix      <- strsplit(textFile[6], " ", fixed = TRUE)[[1]][11]
  gapOpening         <- as.numeric(strsplit(textFile[7], " ", fixed = TRUE)[[1]][14])
  gapExtension       <- as.numeric(strsplit(textFile[8], " ", fixed = TRUE)[[1]][12])
  
  return(c(totalProcessors, skipBadNucleotides, scoringMatrix, gapOpening, gapExtension))

}

# Get the constant for a particular alnalysis run
# It will return in this order: 
# totalProcessors, skipBadNucleotides, scoringMatrix, gapOpening, gapExtension
getAnalysisConstant <- function(AnalysisFolder){
  
  textFile <- readLines(paste(AnalysisFolder,"/AnalysisLog.txt", collapse='', sep=''),encoding="UTF-8")
  
  allInMode         <- as.numeric(strsplit(textFile[2], " ", fixed = TRUE)[[1]][16])
  alignmentDistance <- as.numeric(strsplit(textFile[3], " ", fixed = TRUE)[[1]][8])
  sameCuts          <- as.logical(strsplit(textFile[4], " ", fixed = TRUE)[[1]][17])
  primerDimer       <- as.logical(strsplit(textFile[5], " ", fixed = TRUE)[[1]][14])
  primerConstant    <- as.numeric(strsplit(textFile[6], " ", fixed = TRUE)[[1]][6])
  
  return(c(allInMode, alignmentDistance, sameCuts, primerDimer, primerConstant))
  
}

# Make a summary of the variables used during this combination of runs
makeStatisticsHTML2 <- function (htmlPath, PLOTFOLDER, ANALYSISFOLDER, ALIGNMENTFOLDER, RESULTSFOLDER, SERVERURL){
  
  print(SERVERURL)
  
  # Get the alignment constants
  alignmentVariables <- getAlignmentConstant(RESULTSFOLDER)
  
  # Get the analysis constants
  analysisVariables  <- getAnalysisConstant(ANALYSISFOLDER)
  
  # Get the barcodeDF
  barcodeDF  <- read.table(paste(ALIGNMENTFOLDER, "/barcodeFile_results.txt", sep=''), header=TRUE, sep="\t")
  
  # Get the analysisDF
  analysisDF <- read.table(paste(ANALYSISFOLDER, "/analysis_results.txt", sep=''), header=TRUE, sep="\t")
  
  # Find out the variables
  totalBarcode <- nrow(barcodeDF)
  totalConfig  <- nrow(analysisDF)
  skipBadNucleotides <- alignmentVariables[2]
  scoringMatrix      <- alignmentVariables[3]
  gapOpening         <- alignmentVariables[4]
  gapExtension       <- alignmentVariables[5]
  allInMode          <- analysisVariables[1]
  alignmentDistance  <- analysisVariables[2]
  sameCuts           <- analysisVariables[3]
  primerDimer        <- analysisVariables[4]
  dimerbase          <- analysisVariables[5]
  
  # Create the HTML source file. We are going to write here the HTML result for this function.
  # at the end of this function, this will contain the entire html
  sourceFD <-  file(paste(htmlPath,"/statistics.html",collapse = '',sep = '') , open="w")
  
  
  htmlSource <- paste("<!DOCTYPE html>
                      <html>
                      
                      <head>
                        <link rel='stylesheet' type='text/css' href='style.css'>
                      </head>
                      
                      <body>

                      <div id='indexBodyContainer' class='bodyContainer'>

                        <div id='indexHeader' class='header'>
                        
                          <div id='indexLogo' class='logo'>
                        
                            <img src='Logo.png' alt='I am trap in the server room, please send help!'> 
                        
                          </div>
                        
                        </div>

                        <div id='statisticsRest' class='rest'>
                        
                          <div id='statistics' class='map'>
                        
                            <div id='statisticsBackToParent' class='parentTitle boxTitle'>
                              <div id='linkStatisticsTitle' class='boxElement'>
                                <a href='",SERVERURL,"'> Parent Folder </a>
                              </div>
                            </div>
                        
                            <div id='statisticsTitle' class='indexTitle boxTitle'>
                              <div id='linkStatisticsTitle' class='boxElement'>
                                <a href='index.html'> Index </a>
                              </div>
                            </div>

                            <div id='statisticsTitle' class='statsTitle boxTitle'>
                              <div id='linkStatisticsTitle' class='boxElement'>
                                Statistics
                              </div>
                            </div>
                        
                          </div>

                          <div id='statsTableOfContent' class='TOC'>

                            <a name='Statistics'></a>
                            <h2> STATISTICS </h2>
                            <p> We found a total of", totalBarcode," barcodes.</p>
                            <p> We found a total of", totalConfig ," rows in the configuration file.</p>
                            <h3> Filtering Options </h3>
                            ")
  
  writeLines(htmlSource, sourceFD)
  
  if(skipBadNucleotides==TRUE){
    htmlSource <- "<p> We are skipping sequences with N nucleotides.</p>"
    writeLines(htmlSource, sourceFD)
  }
  
  htmlSource <- paste("
                      
                      <h3> Alignment Options </h3>
                      <p> We are using the scoring matrix: ", scoringMatrix , ".</p>
                      <p> The gap opening penalty is: " , gapOpening , ".</p>
                      <p> The gap extension penalty is: " , gapExtension , ".</p>
                      
                      <h3> Criterias Options </h3>
                      ")
  writeLines(htmlSource, sourceFD)
  
  
  
  if(allInMode == 5){
    htmlSource <- "<p>We don't care if the indels are near or far away from the alignment position.</p>"
    writeLines(htmlSource, sourceFD)
  }
  else{
    if(allInMode == 1){
      htmlSource <- "<p> All indels must be within range of at least one alignment position.</p>"
      writeLines(htmlSource, sourceFD)
    }
    else{
      htmlSource <- "<p> At least one indel must be within range of at least one alignment position.</p>"
      writeLines(htmlSource, sourceFD)
    }
    htmlSource <- paste("<p> The distance from the alignment position is " , alignmentDistance , " nucleotides.</p>")
    writeLines(htmlSource, sourceFD)
  }
  
  
  if(sameCuts == TRUE){
    htmlSource <- "<p>The indels in the forward and reverse reads must be in the same relative position with respect the amplicon.</p>"
    writeLines(htmlSource, sourceFD)
  }
  
  if(primerDimer == TRUE){
    htmlSource <- paste("<p>We are taking care of suspiciously looking primer dimers. Indels that are bigger than the amplicon length - (primers length +",dimerbase,") are considered a primer dimer situation</p>")
    writeLines(htmlSource, sourceFD)
  }
  
  htmlSource <- "<p> <a href='index.html'> <span class='backToIndexText'>Back to Index</span> </a> </p>"
  writeLines(htmlSource, sourceFD)
  
  # End the BODY and HTML tags
  htmlSource <- "</div> </div> </div> </body> </html>"
  writeLines(htmlSource, sourceFD)
  close(sourceFD)
  
}

# Make the barcode HTML page with all the barcode plots
makeBarcodeHTML2 <- function(htmlPath, ANALYSISFOLDER, PLOTFOLDER, barcodeSubDF, plotsBarcodesURL, dataURL, SERVERURL, LINKBY){
  
  # Create the HTML source file. We are going to write here the HTML result for this function.
  # at the end of this function, this will contain the entire html
  sourceFD <-  file(paste(htmlPath,"/barcodes.html",collapse = '',sep = '') , open="w")
  
  # Make the DEFINITIONS part
  htmlSource <- paste("<!DOCTYPE html>
                      <html>
                      
                        <head>
                          <link rel='stylesheet' type='text/css' href='style.css'>
                          <script type='text/javascript' src='scripts.js'></script>
                        </head>
                      
                        <body>

                          <div class='bodyContainer'>
                        
                            <div class='header'>
                        
                              <div class='logo'>
                        
                                <img src='Logo.png' alt='I am trap in the server room, please send help!'> 
                        
                              </div>
                        
                            </div>                    


                            <div class='rest'>
    
                              <div id='barcodesMap' class='map'>
    
                                <div id='barcodesBackToParent' class='parentTitle boxTitle'>
                                  <div id='linkParentTitle' class='boxElement'>
                                    <a href='",SERVERURL,"'> Parent Folder </a>
                                  </div>
                                </div>
    
                                <div id='barcodesTitle' class='indexTitle boxTitle'>
                                  <div id='linkIdexTitle' class='boxElement'>
                                    <a href='index.html'> Index </a>
                                  </div>
                                </div>
    
                                <div id='barcodesTitle' class='barcodesTitle boxTitle'>
                                  <div id='linkBarcodesTitle' class='boxElement'>
                                    Barcodes
                                  </div>
                                </div>

                              </div>

                              <div id='barcodesContent' class='TOC'>

                                <a name='Barcodes'></a>
                                <h2> BARCODES </h2>

                                <div>
                                  <a name='BarcodesDefinitions'></a>
                                  <h3> Definitions </h3>
                                  <p> We have the following plots: </p>
                                  <ul>
                                  <li> <b>Assigned:</b> We have three variables, how many assigned (grey) how many weren't reads (blue) and how many were reads (green). </li>
                                  <li> <b>Filter:</b>  In here we summarize how many sequences we lost to the filter options. If they have N nucleotides they didn't pass (red), otherwise the did pass (green).</li>
                                  </ul>
                                  
                                  <p> We have the following columns: </p>
                                  <ul>
                                  <li> <b>Barcode:   </b> An String ID representing the barcode. </li>
                                  <li> <b>Unnasigned:</b> How many sequences did not have forward or reverse primer. </li>
                                  <li> <b>No cut:    </b> How many sequences had the forward and reverse primer, but didn't pass out cut criteria. </li>
                                  <li> <b>Cut Read:  </b> Sequences with both primers, which passed our cut criteria. </li>
                                  <li> <b>Failed:    </b> How many sequences had a nucleotide label as 'N' in the forward or reverse reads files.</li>
                                  <li> <b>Pass:      </b> How many sequences passed all our filters.  </li>
                                  </ul>
                                </div>
  
                                <div class='TOCContent'>
                                <a name='BarcodeIndex'></a>
                                  <h3> Index </h3>
                                    <ul>")
  
  writeLines(htmlSource, sourceFD)
  
  # We need to adjust the names with spaces so the index looks pretty, lets find maximum of characters we need
  maxLengthName <- 0
  for (i in 1:length(barcodeSubDF)){
    
    if(nchar(as.character(barcodeSubDF[[i]][1,]$"Barcode")) > maxLengthName){
      maxLengthName <- nchar(as.character(barcodeSubDF[[i]][1,]$"Barcode"))
    }
    
  }
  
  # Make the INDEX part
  for (i in 1:length(barcodeSubDF)){
    #start            <- as.character(barcodeSubDF[[i]][1,]$"Barcode")
    #end              <- as.character(barcodeSubDF[[i]][nrow(barcodeSubDF[[i]]), ]$"Barcode")
    
    startLength <- nchar(as.character(barcodeSubDF[[i]][1,]$"Barcode"))
    extraString <- paste(rep("_",(maxLengthName - startLength)),collapse='',sep='')
    
    start <- paste(as.character(barcodeSubDF[[i]][1,]$"Barcode"),extraString,collapse='',sep='')
    end   <- as.character(barcodeSubDF[[i]][nrow(barcodeSubDF[[i]]), ]$"Barcode")
    
    htmlSource <-paste("<li class='barcodeLinkLevel_3'>
                          <a href= '#barcode_",start,"'>
                           <span class='barcodeText'>", start ," ----> ", end ,"</span>
                         </a>
                       </li>",collapse='',sep='')
    
    writeLines(htmlSource, sourceFD)
    
  }
  
  htmlSource <- paste("</ul> </div>")
  writeLines(htmlSource, sourceFD)
  
  
  # Make the TABLE part
  barcodeTableDataPath <- paste(dataURL,         "/barcodeFile_results.txt", collapse = '', sep = '')
  
  htmlSource <-paste("
                     <a name='BarcodesTable'></a>
                     <h3> Barcodes Table </h3>
                     <a href='",barcodeTableDataPath,"'> Barcode Table. </a>
                     ")
  
  writeLines(htmlSource, sourceFD)
  
  
  htmlSource <- paste("       <a name='BarcodesAssignedPlot'></a>
                                  <h3> Barcodes Plot </h3>
                            <div id='barcodesCollection' class='barcodeCollection'>")
  
  writeLines(htmlSource, sourceFD)
  
  for (i in 1:length(barcodeSubDF)){
  
    startLength <- nchar(as.character(barcodeSubDF[[i]][1,]$"Barcode"))
    extraString <- paste(rep("_",(maxLengthName - startLength)),collapse='',sep='')
    
    start <- paste(as.character(barcodeSubDF[[i]][1,]$"Barcode"),extraString,collapse='',sep='')
    end   <- as.character(barcodeSubDF[[i]][nrow(barcodeSubDF[[i]]), ]$"Barcode")
    
    hideButton <- paste("buttonHide_", start, collapse='', sep='')
    showButton <- paste("buttonShow_", start, collapse='', sep='')
    
    assignedPlotPath <- paste(plotsBarcodesURL,"/barcode_assigned_",i,".png", collapse = '', sep = '')
    filterPlotPath   <- paste(plotsBarcodesURL,"/barcode_filters_",i,".png" , collapse = '', sep = '')
    
    assignedDataPath <- paste(dataURL,         "/barcode_assigned_",i,".txt", collapse = '', sep = '')
    filterDataPath   <- paste(dataURL,         "/barcode_filters_",i,".txt" , collapse = '', sep = '')
    
    htmlSource <- paste("

                        <a name= 'barcode_",start,"'></a>
                        <div  class='boxContainer'>
                          <div  class='boxContainerHeader barcodeHeader'>
                            <div class='boxContainerHeaderStart'>
                              ", start, "
                            </div> 
  
                            <div class='boxContainerHeaderTo'>
                              ->
                            </div> 
  
                            <div class='boxContainerHeaderEnd'>
                              ", end, "
                            </div> 
  
                            <div class='boxContainerHeaderVoid'>
                              
                            </div> 
  
                            <div id = '", hideButton,"' class='boxContainerHeaderExpand' onclick='hide(",start,")'>
                              <div> - </div>
                            </div> 
                          </div>

                          <div id = '", start ,"' class='boxContainerBody2'>

                            <div class='boxPlotList2 barcodesPlots'>

                              <div class='boxPlotBox2'>
  
                                <a class='hiddenLink' href=", assignedDataPath, ">
                                <div class='plotTitle'>
                                  <div> Assigned </div> 
                                </div>
                                </a>
                            
                                <div class='plotImage'>
                                   <a href="  , assignedPlotPath,">
                                   <img src=" , assignedPlotPath, " alt='Top 15 no cuts'>
                                   </a>
                                </div>

                              </div> <!-- Close boxPlot -->

                              <div class='boxPlotBox2'>
  
                                <a class='hiddenLink' href=", filterDataPath, ">
                                <div class='plotTitle'>
                                  <div> Filter </div> 
                                </div>
                                </a>
                            
                                <div class='plotImage'>
                                   <a href="  , filterPlotPath,">
                                   <img src=" , filterPlotPath, " alt='Top 15 no cuts'>
                                   </a>
                                </div>

                              </div><!-- Close boxPlot -->


                            </div> <!-- Close the First list with the barcodes plots, assigned and filter -->", sep='', collapse='')
    writeLines(htmlSource, sourceFD)

    
    
    # Lets make the unnasigned div first
    # We are going to make a list of plots. How many list depends on how many barcodes per group we have.
    # We are always going to get 5 images per row.    
    numberOfBarcodes <- nrow(barcodeSubDF[[i]])
    numberOfList <- ceiling (numberOfBarcodes / LINKBY)
    remainder    <- numberOfBarcodes %% LINKBY
    howManyIDs <- 1

    # For each of the list , print a list of links inside that list
    for (j in 1:numberOfList){
  
      htmlSource <- "<div class='boxPlotList2 histogramPlots'>"
      writeLines(htmlSource, sourceFD)
      
      boundary <- 0
            
      # If this is the last group, just write the last links, which is the modulus
      # If this is the last group and the remainder is 0, they match perfectly, so keep the LINKBY
      if(j == numberOfList && remainder !=0){
        boundary <- remainder
      }
      else{
        boundary <- LINKBY
      }  
      
      for (k in 1:boundary){
              
        barcodeID <- as.character(barcodeSubDF[[i]][howManyIDs,]$"Barcode")
              
        histogramLink <- paste(plotsBarcodesURL,"/", barcodeID, "_unnasigned.png", collapse = '', sep = '')
        histogramData <- paste(dataURL,         "/", barcodeID, "_unnasigned.txt", collapse = '', sep = '')
        
        htmlSource <- paste("<div class='boxPlotBox2'>",
                            " <a class='hiddenLink' href=",histogramData,">",
                            "     <div class='plotTitle'>",
                            "         <div>", barcodeID , "</div>",
                            "     </div>
                              </a>
                                     
                              <div class='plotImage'>",
                            "     <a href=",histogramLink,">",
                            "         <img src=",histogramLink," alt='Top 15 unnasigned'>",
                            "     </a>
                              </div>
                            </div>")
        writeLines(htmlSource, sourceFD)
              
        howManyIDs <- howManyIDs + 1
              
      }
      
      
      
      htmlSource <- "</div>"
      writeLines(htmlSource, sourceFD)
        


        
      }
                        
    htmlSource <- "</div> <!-- Close Body -->
                   </div> <!-- Close Container -->"
    writeLines(htmlSource, sourceFD)
    
  }
  
  htmlSource <- paste("       </div> <!-- Close Collection of Containers -->")
  writeLines(htmlSource, sourceFD)
  

  
  htmlSource <- "<p> <a href='index.html'> <span class='backToIndexText'>Back to Index</span> </a> </p>"
  writeLines(htmlSource, sourceFD)
  
  # End the BODY and HTML tags
  htmlSource <- "</div> </div> </div> </body> </html>"
  writeLines(htmlSource, sourceFD)
  close(sourceFD)
  
}

# Make the cutrate HTML page with all the cutrate plots
makeCutRatesHTML3 <- function(htmlPath, ANALYSISFOLDER, PLOTFOLDER, LINKBY, analyisisSubDF,
                              plotsCutratesURL, plotsFrameshiftURL, dataURL, SERVERURL){
  
  # Create the HTML source file. We are going to write here the HTML result for this function.
  # at the end of this function, this will contain the entire html
  sourceFD <-  file(paste(htmlPath,"/frameshift.html",collapse = '',sep = '') , open="w")
  
  htmlSource <- paste("<!DOCTYPE html>
                      <html>

                      <head>
                          <link rel='stylesheet' type='text/css' href='style.css'>
                          <script type='text/javascript' src='scripts.js'></script>
                      </head>

                      <body>
                      
                      <div id='cutratesBodyContainer' class='bodyContainer'>
                      
                        <div id='cutratesHeader' class='header'>
                      
                          <div id='cutratesLogo' class='logo'>
                      
                            <img src='Logo.png' alt='I am trap in the server room, please send help!'> 
                      
                          </div>
                      
                        </div>     
                      
                      <div id='cutratesRest' class='rest'>
                      
                      <div id='cutratesMap' class='map'>
                      
                      <div id='cutratesBackToParent' class='parentTitle boxTitle'>
                      <div id='linkParentTitle' class='boxElement'>
                      <a href='",SERVERURL,"'> Parent Folder </a>
                      </div>
                      </div>
                      
                      <div id='cutratesTitle' class='indexTitle boxTitle'>
                      <div id='linkIdexTitle' class='boxElement'>
                      <a href='index.html'> Index </a>
                      </div>
                      </div>
                      
                      <div  class='analysisTitle boxTitle'>
                        <div  class='boxElement'>
                          Analysis
                        </div>
                      </div>
                      
                      </div>
                      
                      
                      <div id='barcodesContent' class='TOC'>
                      
                      <a name='Cutrates'></a>
                      <h2> CUT RATES </h2>
                      
                      <div>
                      <a name='CutratesDefinitions'></a>
                      <h3> Definitions </h3>
                      <p> We have the following plots: </p>
                      <ul>
                      <li> <b>Normalized:</b> The ratio of sequences that have cuts divided by the number of sequences. </li>
                      <li> <b>Absolute:  </b> How many sequences we have (green + red), and how many sequences are cut reads (green). </li>
                      <li> <b>Frameshift:</b> Wether the deletions we have in the sequence preserve the frame. </li>
                      </ul>
                      
                      <p> We have the following columns: </p>
                      <ul>
                      <li> <b>Cut Reads: </b> These are sequences that have a cut read, hence passed all our criterias. </li>
                      <li> <b>No cut:    </b> These are the sequences that didn't have a cut read. </li>
                      <li> <b>No deletion:   </b> The ratio of sequences that had no deletions. </li>
                      <li> <b>Preserve Frame:</b> The ratio of sequences that have only one deletion of length divisible by 3.</li>
                      <li> <b>Frameshift:    </b> The ratio of sequences that have only one deletion of length not divisible by 3. </li>
                      <li> <b>Multiple:      </b> The ratio of sequences that have more than one deletion. </li>
                      </ul>
                      </div>
                      
                      <div class='TOCContent'>
                      <a name='CutrateIndex'></a>
                      <h3> Index </h3>
                      <ul>")
  writeLines(htmlSource, sourceFD)
  
  # We need to adjust the names with spaces so the index looks pretty, lets find maximum of characters we need
  maxLengthName <- 0
  for (i in 1:length(analyisisSubDF)){
    
    if(nchar(as.character(analyisisSubDF[[i]][1,]$"ID")) > maxLengthName){
      maxLengthName <- nchar(as.character(analyisisSubDF[[i]][1,]$"ID"))
    }
    
  }
  
  for (i in 1:length(analyisisSubDF)){
    
    startLength <- nchar(as.character(analyisisSubDF[[i]][1,]$"ID"))
    extraString <- paste(rep("_",(maxLengthName - startLength)),collapse='',sep='')
    
    start <- paste(as.character(analyisisSubDF[[i]][1,]$"ID"),extraString,collapse='',sep='')
    end   <- as.character(analyisisSubDF[[i]][nrow(analyisisSubDF[[i]]), ]$"ID")
    
    htmlSource <-paste("<li class='cutrateLinkLevel_3'>
                       <a href= '#cutrates_",start,"'>
                       <span class='cutrateText'>", start ," ----> ", end ,"</span>
                       </a>
                       </li>",collapse='',sep='')
    
    writeLines(htmlSource, sourceFD)
    
  }
  
  htmlSource <- paste("</ul> </div>")
  writeLines(htmlSource, sourceFD)
  
  
  htmlSource <- paste("       <a name='cutratePlot'></a>
                      <h3> Plots </h3>
                      <div id='cutratesCollection' class='boxCollection'>")
  
  writeLines(htmlSource, sourceFD)
  
  
  for (i in 1:length(analyisisSubDF)){
    
    startLength <- nchar(as.character(analyisisSubDF[[i]][1,]$"ID"))
    extraString <- paste(rep("_",(maxLengthName - startLength)),collapse='',sep='')
    
    start <- paste(as.character(analyisisSubDF[[i]][1,]$"ID"),extraString,collapse='',sep='')
    end   <- as.character(analyisisSubDF[[i]][nrow(analyisisSubDF[[i]]), ]$"ID")
    
    hideButton <- paste("buttonHide_", start, collapse='', sep='')
    showButton <- paste("buttonShow_", start, collapse='', sep='')
    
    normalPlotPath     <- paste(plotsCutratesURL,  "/Cut_Rates_",i,".png"            , collapse = '', sep = '')
    absolutePlotPath   <- paste(plotsCutratesURL,  "/Cut_Rates_Absolute_",i,".png"   , collapse = '', sep = '')
    frameshiftPlotPath <- paste(plotsFrameshiftURL,"/Frameshift_",i,".png"           , collapse = '', sep = '')
    frameCutsPlotPath  <- paste(plotsFrameshiftURL,"/Frameshift_ReadsOnly_",i,".png" , collapse = '', sep = '')
    
    normalDataPath     <- paste(dataURL,           "/Cut_Rates_",i,".txt"            , collapse = '', sep = '')
    absoluteDataPath   <- paste(dataURL,           "/Cut_Rates_Absolute_",i,".txt"   , collapse = '', sep = '')
    frameshiftDataPath <- paste(dataURL,           "/Frameshift_",i,".txt"           , collapse = '', sep = '')
    frameCutsDataPath  <- paste(dataURL,           "/Frameshift_ReadsOnly_",i,".txt" , collapse = '', sep = '')
    
    htmlSource <- paste("
                        <div id='cutratesContainer' class='boxContainer'>
                          <div id='cutrates_",start,"' class='boxContainerHeader experimentHeader'>
                            <div id='cutratesContainerHeaderStart' class='boxContainerHeaderStart'>
                              ", start, "
                            </div> 
                        
                            <div id='cutratesContainerHeaderTo' class='boxContainerHeaderTo'>
                              ->
                            </div> 
                        
                            <div id='cutratesContainerHeaderEnd' class='boxContainerHeaderEnd'>
                              ", end, "
                            </div> 
                        
                            <div id='cutratesContainerHeaderVoid' class='boxContainerHeaderVoid'>
                        
                            </div> 
  
                            <div id = '", hideButton,"' class='boxContainerHeaderExpand' onclick='hide(",start,")'>
                              <div> - </div>
                            </div> 

                          </div>
                        
                          <div id= '", start ,"' class='boxContainerBody'>
                        
                            <div class='boxPlotList'>

                              <div class='boxPlotBox'>

                                <a class='hiddenLink' href=", normalDataPath, ">
                                <div class='plotTitle'>
                                  <div> Normalize </div>
                                </div>
                                </a>
                          
                                <div class='plotImage'>
                                  <a href="  , normalPlotPath,">
                                  <img src=" , normalPlotPath, " alt='Assigned sequences for this barcodes' >
                                  </a>
                                </div>

                              </div>

                              <div class='boxPlotBox'>

                                <a class='hiddenLink' href=", frameshiftDataPath, ">
                                <div class='plotTitle'>
                                  <div> Frameshift </div> 
                                </div>
                                </a>
                          
                                <div class='plotImage'>
                                  <a href="  , frameshiftPlotPath,">
                                  <img src=" , frameshiftPlotPath, " alt='Assigned sequences for this barcodes' >
                                  </a>
                                </div>

                              </div>

                            </div>

                            <div class='boxPlotList'>

                              <div class='boxPlotBox'>

                                <a class='hiddenLink' href=", absoluteDataPath, ">
                                <div class='plotTitle'>
                                  <div> Absolute </div> 
                                </div>
                                </a>
                          
                                <div class='plotImage'>
                                  <a href="  , absolutePlotPath, ">
                                  <img src=" , absolutePlotPath, " alt='Cut Rate, absolute numbers' >
                                  </a>
                                </div>

                              </div>

                              <div>

                                <a class='hiddenLink' href=", frameCutsDataPath, ">
                                <div class='plotTitle'>
                                  <div> Frameshift (Cuts Only) </div> 
                                </div>
                                </a>
                          
                                <div class='plotImage'>
                                  <a href="  , frameCutsPlotPath,">
                                  <img src=" , frameCutsPlotPath, " alt='Frameshift, only reads' >
                                  </a>
                                </div>

                              </div>

                            </div>

                          </div>

                          <div id= 'LinkCollection' class='linkBox'>

                        ", sep='', collapse='')
    writeLines(htmlSource, sourceFD)

    # Now we are going to make the link box at the end of the container
    # We just need to divide the GROUPBY number of links into groups of size LINKBY
    
    numberOfIDs <- nrow(analyisisSubDF[[i]])
    
    numberOfList <- ceiling (numberOfIDs / LINKBY)
    remainder    <- numberOfIDs %% LINKBY
    howManyIDs <- 1
    
    # For each of the list , print a list of links inside that list
    for (j in 1:numberOfList){
      
      boundary <- 0
      
      # If this is the last group, just write the last links, which is the modulus
      # If this is the last group and the remainder is 0, they match perfectly, so keep the LINKBY
      if(j == numberOfList && remainder !=0){
        boundary <- remainder
      }
      else{
        boundary <- LINKBY
      }
      
      htmlSource <-  "<div class='linkList'>"
      writeLines(htmlSource, sourceFD)
      
      for (k in 1:boundary){
      
        geneID <- as.character(analyisisSubDF[[i]][howManyIDs,]$"ID")

        
        htmlSource <- paste("
                            <div class='linkContainer'>
                            
                            <a href= 'htmlResults_",i,".html#",geneID, "'>", geneID , "</a>

                            </div>
                            ", sep='', collapse='')
        writeLines(htmlSource, sourceFD)
        
        howManyIDs <- howManyIDs + 1
        
      }

      htmlSource <-  "</div>"
      writeLines(htmlSource, sourceFD)
      
    }
    

    # Close for the link box and the whole container
    htmlSource <- paste("

                          </div>
                        
                        </div> 
                        ", sep='', collapse='')
    writeLines(htmlSource, sourceFD)
    
  }
  
  
  htmlSource <- paste("       </div>")
  writeLines(htmlSource, sourceFD)
  
  
  htmlSource <-paste("
                     <a name='Cut Rate Table'></a>
                     <h3> Analysis info (TAR.GZ) </h3>
                     <a href='http://www.ii.uib.no/~eivindv/Jason2/Analysis.tar.gz'> Barcode Table. </a>
                     ")
  
  writeLines(htmlSource, sourceFD)
  
  htmlSource <- "<p> <a href='index.html'> <span class='backToIndexText'>Back to Index</span> </a> </p>"
  writeLines(htmlSource, sourceFD)
  
  # End the BODY and HTML tags
  htmlSource <- "</div> </div> </div> </body> </html>"
  writeLines(htmlSource, sourceFD)
  close(sourceFD)
  
}

# Make the Resutls HTML page with all the cutrate plots
makeResultsHTML3 <- function(htmlPath, ANALYSISFOLDER, PLOTFOLDER, analyisisSubDF,
                             PLOTTARGET, dataURL){
  
  # We need to make one results for each group
  for (i in 1:length(analyisisSubDF)){
    
    start            <- as.character(analyisisSubDF[[i]][1,]$"ID")
    end              <- as.character(analyisisSubDF[[i]][nrow(analyisisSubDF[[i]]), ]$"ID")
    
    # Create the HTML source file. We are going to write here the HTML result for this function.
    # at the end of this function, this will contain the entire html
    sourceFD <-  file(paste(htmlPath,"/htmlResults_",i,".html",collapse = '',sep = '') , open="w")
    
    htmlSource <- paste("<!DOCTYPE html>
                        <html>
                        
                        <head>
                        <link rel='stylesheet' type='text/css' href='style.css'>
                        </head>
                        
                        <body>
                        
                        <div id='cutratesBodyContainer' class='bodyContainer'>
                        
                        <div id='cutratesHeader' class='header'>
                        
                        <div id='cutratesLogo' class='logo'>
                        
                        <img src='Logo.png' alt='I am trap in the server room, please send help!'> 
                        
                        </div>
                        
                        </div> 
                        
                        <div class='rest'>
                        
                        <div class='map'>
                        
                        <div class='parentTitle boxTitle'>
                        <div id='linkParentTitle' class='boxElement'>
                        <a href='",SERVERURL,"'> Parent Folder </a>
                        </div>
                        </div>
                        
                        <div class='indexTitle boxTitle'>
                        <div class='boxElement'>
                        <a href='index.html'> Index </a>
                        </div>
                        </div>
                        
                        <div  class='resultsTitle boxTitle'>
                        <div  class='boxElement'>
                        Results
                        </div>
                        </div>
                        
                        </div>
                        
                        <div id='barcodesContent' class='TOC'>
                        
                        <a name='Results'></a>
                        <h2> Results from ", start, " to ", end, " </h2>
                        
                        <div>
                        <a name='Definitions'></a>
                        <h3> Definitions </h3>
                        <p> We have the following plots: </p>
                        <ul>
                        <li> <b>Archplot:        </b> Plot that shows the deletions and cuts respect the amplicon. </li>
                        <li> <b>Mutations:       </b> Plot that shows the mutations for each nucleotide in the amplicon. </li>
                        <li> <b>Top15 Cuts:      </b> Plot that shows the most abundant reads with cuts for that particular Gene ID. </li>
                        <li> <b>Top15 Sequences: </b> Plot that shows the most abundant reads (with or without cuts) for that particular Gene ID. </li>
                        </ul>
                        
                        <p> We have the following columns: </p>
                        <ul>
                        <li> <b>Deletion: </b> A deletion in the alignment from the amplicon coordinates. </li>
                        <li> <b>Cut:      </b> A cut (deletion that passed the criteria) from the amplicon coordinates. </li>
                        <li> <b>Frequency:</b> The number of sequences divided by the total of sequences that have that particular cut or deletion. </li>
                        </ul>
                        </div>
                        
                        <div class='TOCContent'>
                        <a name='Index'></a>
                        <h3> Index </h3>
                        <ul>")
    
    writeLines(htmlSource, sourceFD)
    
    # Make the index for each element of the table
    for (j in 1:nrow(analyisisSubDF[[i]])){
      
      id      <- analyisisSubDF[[i]][j,]$"ID"
      
      htmlSource <- paste("
                          <li> <a href=#", id , ">", id, "</a> </li>
                          ", sep='', collapse='')
      writeLines(htmlSource, sourceFD)
      
    }
    
    htmlSource <- paste("</ul> </div>")
    writeLines(htmlSource, sourceFD)
    
    
    htmlSource <- paste("       <a name='cutratePlot'></a>
                        <h3> Plots </h3>
                        <div id='cutratesCollection' class='boxCollection'>")
    
    writeLines(htmlSource, sourceFD)
    
    # Make each of the Box Containers
    for (j in 1:nrow(analyisisSubDF[[i]])){
      
      id      <- analyisisSubDF[[i]][j,]$"ID"
      barcode <- analyisisSubDF[[i]][j,]$"Barcode"
      
      hideButton <- paste("buttonHide_", start, collapse='', sep='')
      showButton <- paste("buttonShow_", start, collapse='', sep='')
      
      archplotName   <- paste(id,"_",barcode,"_archplot.png", collapse = '', sep = '')
      histrogramName <- paste(id,"_",barcode,"_histogramBOTHCUT.png", collapse = '', sep = '')
      
      archPlotPath       <- paste(PLOTFOLDER,"/", archplotName , collapse = '', sep = '')      
      histrogramPlotPath <- paste(PLOTFOLDER,"/", histrogramName , collapse = '', sep = '')
      
      archPlotTarget       <- paste(PLOTTARGET,"/", archplotName , collapse = '', sep = '')      
      histrogramPlotTarget <- paste(PLOTTARGET,"/", histrogramName , collapse = '', sep = '')      
      
      #archPlotTDString  <- ""
      #histogramTDString <- ""
      
      # Get the images with links that goes inside the plot containers.
      # They might not exist, we generate a text message in that case
      if(file.exists(archPlotPath)){
        archPlotTDString <- paste("
                                  <a href="  , archPlotTarget,">
                                  <img src=" , archPlotTarget, " alt='Archplot'>
                                  </a>")
      }
      else{
        archPlotTDString <- "Archplot not generated"
      }
      
      if(file.exists(histrogramPlotPath)){
        histogramTDString <- paste("
                                   <a href="  , histrogramPlotTarget,">
                                   <img src=" , histrogramPlotTarget, " alt='Top 15 with cuts'>
                                   </a>")
      }
      else{
        histogramTDString <- "<div> Histogram not generated </div>"
      }
      
      htmlSource <- paste("
                          <div id='", id," class='boxContainer'>
                          
                          <div class='boxContainerHeader resultsHeader'>
                            <div class='boxContainerHeaderStart'>
                              ", id, "
                            </div> 
                          
                            <div  class='boxContainerHeaderTo'>
                              |
                            </div> 
                          
                            <div  class='boxContainerHeaderEnd'>
                              ", barcode, "
                            </div> 
                          
                            <div  class='boxContainerHeaderLink'>
                              | <a href='#Index'> Index </a>
                            </div> 
                          
                            <div  class='boxContainerHeaderVoid'>
                          
                            </div> 
                          
                            <div id = '", hideButton,"' class='boxContainerHeaderExpand' onclick='hide(",start,")'>
                              <div> - </div>
                            </div> 
                          
                          </div>
                          
                          <div  id = '",start,"' class='boxContainerBody'>
                          
                            <div class='boxPlotList ampliconPlots'>
                          
                              <div class='boxPlotBox'>
                                <div class='plotTitle'>
                                  <div> Archplot </div> 
                                </div>

                                <div class='plotImage'>
                                  ",archPlotTDString,"
                                </div>
                              </div>"
                          
#                               <div class='boxPlotBox'>
#                                 <div class='plotTitle'>
#                                   <div> Mutations </div> 
#                                 </div>
#                           
#                                 <div class='plotImage'>
#                                   <p>The ghost of plots yet to come.</p>
#                                 </div>
#                               </div>
                          
                          ,"</div>
                          
                          <div class='boxPlotList histogramPlots'>
                          
                            <div class='boxPlotBox'>
                              <div class='plotTitle'>
                                <div> Top15 Reads </div> 
                              </div>
                          
                              <div class='plotImage'>
                                ",histogramTDString,"
                              </div>
                            </div>"
                          
#                             <div class='boxPlotBox'>
#                               <div class='plotTitle'>
#                                 <div> Top 15 Sequences </div> 
#                               </div>
#                           
#                               <div class='plotImage'>
#                                 <p>The ghost of plots yet to come.</p>
#                               </div>
#                             </div>
                          
                          ,"</div> <!--Close the Plot List-->
                          
                          </div> <!--Close the Body Container-->
                          
                          </div> <!--Close the Container-->
                          
                          ", sep='', collapse='')
      writeLines(htmlSource, sourceFD)
      
      
      }
    
    # Close the Box Collection
    htmlSource <- paste("       </div>")
    writeLines(htmlSource, sourceFD)
    
    htmlSource <- "<p> <a href='index.html'> <span class='backToIndexText'>Back to Index</span> </a> </p>"
    writeLines(htmlSource, sourceFD)
    
    # End the BODY and HTML tags
    htmlSource <- "</div> </div> </div> </body> </html>"
    writeLines(htmlSource, sourceFD)
    close(sourceFD)
    
    }
  
}

# Make the Resutls HTML page with all the cutrate plots
makeResultsHTML4 <- function(htmlPath, ANALYSISFOLDER, PLOTFOLDER, analyisisSubDF,
                             PLOTTARGET, dataURL, SERVERURL){
  
  # We need to make one results for each group
  for (i in 1:length(analyisisSubDF)){
    
    start            <- as.character(analyisisSubDF[[i]][1,]$"ID")
    end              <- as.character(analyisisSubDF[[i]][nrow(analyisisSubDF[[i]]), ]$"ID")
    
    # Create the HTML source file. We are going to write here the HTML result for this function.
    # at the end of this function, this will contain the entire html
    sourceFD <-  file(paste(htmlPath,"/htmlResults_",i,".html",collapse = '',sep = '') , open="w")
    
    htmlSource <- paste("<!DOCTYPE html>
                        <html>
                        
                        <head>
                          <link rel='stylesheet' type='text/css' href='style.css'>
                          <script type='text/javascript' src='scripts.js'></script>
                        </head>
                        
                        <body>
                        
                        <div class='bodyContainer'>
                        
                        <div class='header'>
                        
                        <div class='logo'>
                        
                        <img src='Logo.png' alt='I am trap in the server room, please send help!'> 
                        
                        </div>
                        
                        </div> 
                        
                        <div class='rest'>
                        
                        <div class='map'>
                        
                        <div class='parentTitle boxTitle'>
                        <div id='linkParentTitle' class='boxElement'>
                        <a href='",SERVERURL,"'> Parent Folder </a>
                        </div>
                        </div>
                        
                        <div class='indexTitle boxTitle'>
                        <div class='boxElement'>
                        <a href='index.html'> Index </a>
                        </div>
                        </div>
                        
                        <div  class='resultsTitle boxTitle'>
                        <div  class='boxElement'>
                        Results
                        </div>
                        </div>
                        
                        </div>
                        
                        <div id='barcodesContent' class='TOC'>
                        
                        <a name='Results'></a>
                        <h2> Results from ", start, " to ", end, " </h2>
                        
                        <div>
                        <a name='Definitions'></a>
                        <h3> Definitions </h3>
                        <p> We have the following plots: </p>
                        <ul>
                        <li> <b>Archplot:        </b> Plot that shows the deletions and cuts respect the amplicon. </li>
                        <li> <b>Mutations:       </b> Plot that shows the mutations for each nucleotide in the amplicon. </li>
                        <li> <b>Top15 Cuts:      </b> Plot that shows the most abundant reads with cuts for that particular Gene ID. </li>
                        <li> <b>Top15 Sequences: </b> Plot that shows the most abundant reads (with or without cuts) for that particular Gene ID. </li>
                        </ul>
                        
                        <p> We have the following columns: </p>
                        <ul>
                        <li> <b>Deletion: </b> A deletion in the alignment from the amplicon coordinates. </li>
                        <li> <b>Cut:      </b> A cut (deletion that passed the criteria) from the amplicon coordinates. </li>
                        <li> <b>Frequency:</b> The number of sequences divided by the total of sequences that have that particular cut or deletion. </li>
                        </ul>
                        </div>
                        
                        <div class='TOCContent'>
                        <a name='Index'></a>
                        <h3> Index </h3>
                        <ul>")
    
    writeLines(htmlSource, sourceFD)
    
    # Make the index for each element of the table
    for (j in 1:nrow(analyisisSubDF[[i]])){
      
      id      <- analyisisSubDF[[i]][j,]$"ID"
      
      htmlSource <- paste("
                          <li> <a href=#", id , ">", id, "</a> </li>
                          ", sep='', collapse='')
      writeLines(htmlSource, sourceFD)
      
    }
    
    htmlSource <- paste("</ul> </div>")
    writeLines(htmlSource, sourceFD)
    
    
    htmlSource <- paste("       <a name='cutratePlot'></a>
                        <h3> Plots </h3>
                        <div id='cutratesCollection' class='boxCollection'>")
    
    writeLines(htmlSource, sourceFD)
    
    # Make each of the Box Containers
    for (j in 1:nrow(analyisisSubDF[[i]])){
      
      id      <- analyisisSubDF[[i]][j,]$"ID"
      barcode <- analyisisSubDF[[i]][j,]$"Barcode"
      
      idBoxContainer <- paste(id,"_boxContainer",collapse="",sep="")
      
      hideButton <- paste("buttonHide_", start, collapse='', sep='')
      showButton <- paste("buttonShow_", start, collapse='', sep='')
      
      archplotName   <- paste(id,"_",barcode,"_archplot.png", collapse = '', sep = '')
      histrogramName <- paste(id,"_",barcode,"_histogramBOTHCUT.png", collapse = '', sep = '')
      histrogramDELName <- paste(id,"_",barcode,"_histogramNOCUT.png", collapse = '', sep = '')
      histrogramCUTName <- paste(id,"_",barcode,"_histogramGOODCUT.png", collapse = '', sep = '')
      
      archPlotPath       <- paste(PLOTFOLDER,"/", archplotName , collapse = '', sep = '')      
      histrogramPlotPath <- paste(PLOTFOLDER,"/", histrogramName , collapse = '', sep = '')
      histrogramDELPlotPath <- paste(PLOTFOLDER,"/", histrogramDELName , collapse = '', sep = '')
      histrogramCUTPlotPath <- paste(PLOTFOLDER,"/", histrogramCUTName , collapse = '', sep = '')
      
      archPlotTarget       <- paste(PLOTTARGET,"/", archplotName , collapse = '', sep = '')      
      histrogramPlotTarget <- paste(PLOTTARGET,"/", histrogramName , collapse = '', sep = '')      
      histrogramDELPlotTarget <- paste(PLOTTARGET,"/", histrogramDELName , collapse = '', sep = '')      
      histrogramCUTPlotTarget <- paste(PLOTTARGET,"/", histrogramCUTName , collapse = '', sep = '')      
      
      archDataName          <- paste(id,"_",barcode,"_archplot.txt", collapse = '', sep = '')
      histrogramDataName    <- paste(id,"_",barcode,"_histogramBOTHCUT.txt", collapse = '', sep = '')
      histrogramDELDataName <- paste(id,"_",barcode,"_histogramNOCUT.txt", collapse = '', sep = '')
      histrogramCUTDataName <- paste(id,"_",barcode,"_histogramGOODCUT.txt", collapse = '', sep = '')

      archDataTarget          <- paste(dataURL,"/", archDataName , collapse = '', sep = '')      
      histrogramDataTarget    <- paste(dataURL,"/", histrogramDataName , collapse = '', sep = '')      
      histrogramDELDataTarget <- paste(dataURL,"/", histrogramDELDataName , collapse = '', sep = '')      
      histrogramCUTDataTarget <- paste(dataURL,"/", histrogramCUTDataName , collapse = '', sep = '')      
      
      
      archPlotTDString  <- ""
      histogramTDString <- ""
      histogramTDDELString <- ""
      histogramTDCUTString <- ""
      
      # Get the images with links that goes inside the plot containers.
      # They might not exist, we generate a text message in that case
      if(file.exists(archPlotPath)){
        archPlotTDString <- paste("
                                  <a href="  , archPlotTarget,">
                                  <img src=" , archPlotTarget, " alt='Archplot'>
                                  </a>")
      }
      else{
        archPlotTDString <- "Archplot not generated"
      }
      
      if(file.exists(histrogramPlotPath)){
        histogramTDString <- paste("
                                   <a href="  , histrogramPlotTarget,">
                                   <img src=" , histrogramPlotTarget, " alt='Top 15 with cuts'>
                                   </a>")
      }
      else{
        histogramTDString <- "<div> Histogram not generated </div>"
      }
      
      if(file.exists(histrogramDELPlotPath)){
        histogramTDDELString <- paste("
                                   <a href="  , histrogramDELPlotTarget,">
                                   <img src=" , histrogramDELPlotTarget, " alt='Top 15 no cuts'>
                                   </a>")
      }
      else{
        histogramTDDELString <- "<div> Histogram not generated </div>"
      }
      
      if(file.exists(histrogramCUTPlotPath)){
        histogramTDCUTString <- paste("
                                   <a href="  , histrogramCUTPlotTarget,">
                                   <img src=" , histrogramCUTPlotTarget, " alt='Top 15 with cuts'>
                                   </a>")
      }
      else{
        histogramTDCUTString <- "<div> Histogram not generated </div>"
      }
      
      htmlSource <- paste("
                          <a name=",id,"></a>
                          <div id='", id," class='boxContainer'>
                          
                          <div class='boxContainerHeader resultsHeader'>

                            <div class='boxContainerHeaderStart'>
                            ", id, "
                            </div> 
                          
                            <div  class='boxContainerHeaderTo'>
                            |
                            </div> 
                          
                            <div  class='boxContainerHeaderEnd'>
                            ", barcode, "
                            </div> 
                          
                            <div  class='boxContainerHeaderLink'>
                                | <a href='#Index'> Index </a>
                            </div> 
                          
                            <div  class='boxContainerHeaderVoid'>
                            
                            </div> 

                            <div id = '", hideButton,"' class='boxContainerHeaderExpand' onclick='hide(",idBoxContainer,")'>
                                <div> - </div>
                            </div> 
                          
                          </div>
                          
                          
                          <div  id = '",idBoxContainer,"' class='boxContainerBody2'>
                          
                          <div class='boxPlotList2 ampliconPlots'>
                          
                            <div class='boxPlotBox2'>
                              <a class='hiddenLink' href=", archDataTarget, ">
                              <div class='plotTitle'>
                                <div> Archplot </div> 
                              </div>
                              </a>
                          
                              <div class='plotImage'>
                                ",archPlotTDString,"
                              </div>
                            </div>

                          </div>
                          
                          <div class='boxPlotList2 histogramPlots'>
                          
                            <div class='boxPlotBox2'>

                              <a class='hiddenLink' href=", histrogramDataTarget, ">
                              <div class='plotTitle'>
                                <div> Top15 Both </div> 
                              </div>
                              </a>
                          
                              <div class='plotImage'>
                                ",histogramTDString,"
                              </div>
                            </div>

                            <div class='boxPlotBox2'>

                              <a class='hiddenLink' href=", histrogramCUTDataTarget, ">
                              <div class='plotTitle'>
                                <div> Top15 With Cuts </div> 
                              </div>
                              </a>
                          
                              <div class='plotImage'>
                                ",histogramTDCUTString,"
                              </div>
                            </div>

                            <div class='boxPlotBox2'>

                              <a class='hiddenLink' href=", histrogramDELDataTarget, ">
                              <div class='plotTitle'>
                                <div> Top15 No Cuts </div> 
                              </div>
                              </a>
                          
                              <div class='plotImage'>
                                ",histogramTDDELString,"
                              </div>
                            </div>

                          </div> <!--Close the Plot List-->
                          
                          </div> <!--Close the Body Container-->
                          
                          </div> <!--Close the Container-->
                          
                          ", sep='', collapse='')
      writeLines(htmlSource, sourceFD)
      
      
      }
    
    # Close the Box Collection
    htmlSource <- paste("       </div>")
    writeLines(htmlSource, sourceFD)
    
    htmlSource <- "<p> <a href='index.html'> <span class='backToIndexText'>Back to Index</span> </a> </p>"
    writeLines(htmlSource, sourceFD)
    
    # End the BODY and HTML tags
    htmlSource <- "</div> </div> </div> </body> </html>"
    writeLines(htmlSource, sourceFD)
    close(sourceFD)
    
    }
  
}
