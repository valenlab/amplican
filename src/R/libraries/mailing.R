library("sendmailR")
library("png") # HTML functions will use this to resize RAW pngs into tiny previews
library("raster")
#setwd("~/Desktop/git/ampliCan/src/R")

# Constant definitions

# This strings represent part of the names of some folders.
# The program look for this strings in a result folder in order to find each of the subfolders/files
PLOTFOLDERSTRING <- "plots"
RESULTFOLDERSTRING <- "results"
LOGFILESTRING <- "myFinalLog"

# Definitions for each of the columns


sendmail_options(smtpServer="smtp.uib.no")


sendReport <- function (from, to, filePath, htmlPath){
  
  body <- c("Here are the reports",
          mime_part(x=htmlPath,name="config_result.HTML"),
          mime_part(x=filePath,name="config_result.txt")
          )
  
  
  
  to <- "rafanozal@gmail.com"
  subject <- "[ampliCan] Automatic result report"
  
  sendmail(from, to, subject, body)
  

}

dataframeToHTML <- function(dataframeFilePath, htmlBasePath, finalNamePath){
  
  print("HTML")
  
  #Open the dataframe file and put it into a dataframe
  dataframe <- read.table(dataframeFilePath, header=TRUE)
  
  htmlFD <-  file(htmlBasePath, open="r")
  htmlEnding <- "</table> </body> </html>"
  htmlStart <- readLines(htmlBasePath,encoding="UTF-8")
  htmlStart <- paste(htmlStart,collapse='')
  
  columnNames <- colnames(dataframe)
  totalRows <- nrow(dataframe)
  totalColumns <- ncol(dataframe)
  
  # First, lets make the header
  htmlStart <- paste(htmlStart , "<tr>",collapse='')
    
  for(j in 1:totalColumns){
  
    htmlStart <- paste(htmlStart, "<th>",columnNames[j]," </th> ",collapse='')
    
  }
  
  htmlStart <- paste(htmlStart,"</tr>",collapse='')
  
  # Now, lets make each row
  for(i in 1:totalRows){
    htmlStart <- paste(htmlStart , "<tr>",collapse='')
  
    for(j in 1:totalColumns){
      
      htmlStart <- paste(htmlStart, "<td>",toString(dataframe[i,j])," </td>", collapse='')
    }
    
    htmlStart <- paste(htmlStart," </tr>",collapse='')
  }
  
  # Close the html
  htmlStart <- paste(htmlStart, htmlEnding, collapse='')
  
  # Finally; write the results into disk
  resultFD <-  file(finalNamePath, open="w")
  writeLines(c(toString(htmlStart)), resultFD)
  close(resultFD)
  close(htmlFD)
  
  return(htmlStart)

}


generateHTMLReport <- function(folderPath){

  print("Generating Report")



}

resultFolderToHTML <- function(pathToFolder, htmlBasePath){

  #Create a folder where we are going to put every HTML related stuff
  dir.create(file.path(pathToFolder, "/HTML"), showWarnings = FALSE)
  htmlFolder <- paste(pathToFolder,"/HTML",sep = '')
  
  # Prepare the beggining of the web, this is, the HTML headers, the CSS and the Javascript. This part is constant
  # for every results, 
  htmlFD <-  file(htmlBasePath, open="r")
  htmlEnding <- "</body> </html>"
  
  # In here we are going to write the html code.
  # All this function will add text to this variable until the HTML5 document is finish
  htmlSource <- readLines(htmlBasePath,encoding="UTF-8")
  htmlSource <- paste(htmlSource,collapse='')
  
  close(htmlFD) # We don't need the base anymore after this
  
  # Look for the plot and results folder
  plotFolderPath <- list.dirs(pathToFolder, recursive = FALSE)[1]
  resultFolderPath <- list.dirs(pathToFolder, recursive = FALSE)[2]
  
  # Look for the master log
  
  # Load the config result and the barcode results
  configResultPath <- paste(resultFolderPath,"config_results.txt",sep='/')
  barcodeResultPath <- paste(resultFolderPath,"barcode_results.txt",sep='/')
  
  configResultTable <- read.table(configResultPath, header=TRUE, sep="\t")
  barcodeResultTable <- read.table(barcodeResultPath, header=TRUE, sep="\t")
  
  # Make four primary links:
  # --- Link for the statistics results
  # --- Link for the barcodes results
  # ------ Definitions
  # ------ Assigned Plot
  # ------ Filter Plot
  # ------ Barcodes Table
  # --- Link for the config table result
  # ------ Definitions
  # ------ Cut Rate Plot  
  # ------ Frameshift Plot  
  # ------ Config index
  # --- Link for the log text result
  # ------ RAW Log results
  # --- Link for the alignment overview
  
  htmlSource <- paste(htmlSource,"",collapse='')
  
  htmlSource <- paste(htmlSource,"<a name='index'></a> ",collapse='') # Anchor to the index (for back to top buttons)
  htmlSource <- paste(htmlSource,"<div id='pageIndex' class='pageIndex'>",collapse='')
  htmlSource <- paste(htmlSource,"<div id='pageIndexTitle' class='indexTitle'> <h2>Contents</h2> </div>",collapse='')
  htmlSource <- paste(htmlSource,"<div id='pageIndexBody'  class='indexBody'>",collapse='')
  htmlSource <- paste(htmlSource,"<ul>",collapse='')
  
  htmlSource <- paste(htmlSource,'<li class="indexLinkLevel_1"> <a href="#Statistics">   <span class="indexNumber">1</span> <span class="indexText">Statistics</span>   </a> </li>',collapse='')
  
  htmlSource <- paste(htmlSource,'<li class="indexLinkLevel_1"> <a href="#Barcodes">     <span class="indexNumber">2</span> <span class="indexText">Barcodes</span>     </a> 
                                  <ul>
                                  <li class="indexLinkLevel_2"> <a href="#BarcodesDefinitions">   <span class="indexNumber">2</span> <span class="indexText">Definitions</span>   </a></li>
                                  <li class="indexLinkLevel_2"> <a href="#BarcodesAssignedPlot">  <span class="indexNumber">2</span> <span class="indexText">Assigned Plot</span> </a></li>
                                  <li class="indexLinkLevel_2"> <a href="#BarcodesFilterPlot">    <span class="indexNumber">2</span> <span class="indexText">Filter Plot</span>   </a></li>
                                  <li class="indexLinkLevel_2"> <a href="#BarcodesTable">         <span class="indexNumber">2</span> <span class="indexText">Barcode Table</span> </a></li>
                                  </ul>
                                  </li>',collapse='')
  
  htmlSource <- paste(htmlSource,'<li class="indexLinkLevel_1"> <a href="#Config">  <span class="indexNumber">3</span> <span class="indexText">Config Table</span> </a>
                                  <ul>
                                  <li class="indexLinkLevel_2"> <a href="#ConfigDefinitions">    <span class="indexNumber">2</span> <span class="indexText">Definitions</span>   </a></li>
                                  <li class="indexLinkLevel_2"> <a href="#ConfigCutRatePlot">    <span class="indexNumber">2</span> <span class="indexText">Cut Rate Plot</span> </a></li>
                                  <li class="indexLinkLevel_2"> <a href="#ConfigFrameshiftPlot"> <span class="indexNumber">2</span> <span class="indexText">Frameshift Plot</span>   </a></li>
                                  <li class="indexLinkLevel_2"> <a href="#ConfigTable">          <span class="indexNumber">2</span> <span class="indexText">Config Table</span> </a></li>
                                  </ul>
                                  </li>',collapse='')
  
  
  
  
  htmlSource <- paste(htmlSource,'<li class="indexLinkLevel_1"> <a href="#Logs">         <span class="indexNumber">4</span> <span class="indexText">Log</span>          </a> </li>',collapse='')
  htmlSource <- paste(htmlSource,'<li class="indexLinkLevel_1"> <a href="#Alignments">         <span class="indexNumber">5</span> <span class="indexText">Alignments</span>          </a> </li>',collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"",collapse='')
  htmlSource <- paste(htmlSource,"<p></p>",collapse='')
  
  htmlSource <- paste(htmlSource,"</div>",collapse='') # /indexBody
  htmlSource <- paste(htmlSource,"</div>",collapse='') # /index

  # -----------------------------
  # STATISTICS
  # Start the page; give a brief introduction about statistics of the document
  # -----------------------------  
  {
  htmlSource <- paste(htmlSource,"<a name='Statistics'></a> ",collapse='')
  htmlSource <- paste(htmlSource,"<h2> STATISTICS </h2>",collapse='')
  htmlSource <- paste(htmlSource,"<p> We found a total of", nrow(barcodeResultTable)," barcodes.</p>",collapse='')
  htmlSource <- paste(htmlSource,"<p> We found a total of", nrow(configResultTable)," rows in the configuration file.</p>",collapse='')
  htmlSource <- paste(htmlSource,"<h3> Filtering Options </h3>",collapse='')
  if(SKIP_BAD_NUCLEOTIDES==TRUE){
    htmlSource <- paste(htmlSource,"<p> We are skipping sequences with N nucleotides.</p>",collapse='')
  }
  
  htmlSource <- paste(htmlSource,"<h3> Alignment Options </h3>",collapse='')
  htmlSource <- paste(htmlSource,"<p> We are using the scoring matrix: ",SCORING_MATRIX,".</p>",collapse='')
  htmlSource <- paste(htmlSource,"<p> The gap opening penalty is: ",GAP_OPENING,".</p>",collapse='')
  htmlSource <- paste(htmlSource,"<p> The gap extension penalty is: ",GAP_EXTENSION,".</p>",collapse='')
  
  htmlSource <- paste(htmlSource,"<h3> Criterias Options </h3>",collapse='')
  if(ALL_IN_MODE == 5){
    htmlSource <- paste(htmlSource,"<p>We don't care if the indels are near or far away from the alignment position.</p>",collapse='')
  }
  else{
    if(ALL_IN_MODE == 1){
      htmlSource <- paste(htmlSource,"<p> All indels must be within range of at least one alignment position.</p>",collapse='')
    }
    else{
      htmlSource <- paste(htmlSource,"<p> At least one indel must be within range of at least one alignment position.</p>",collapse='')
    }
    htmlSource <- paste(htmlSource,"<p> The distance from the alignment position is ",ALIGNMENT_DISTANCE," nucleotides.</p>",collapse='')
  }
  if(SAME_CUTS == TRUE){
    htmlSource <- paste(htmlSource,"<p>The indels in the forward and reverse reads must be in the same relative position with respect the amplicon.</p>",collapse='')
  }
  if(PRIMERDIMER == TRUE){
    htmlSource <- paste(htmlSource,"<p>We are taking care of suspiciously looking primer dimers. Indels that are bigger than the amplicon length - (primers length +",DIMERBASE," are considered a primer dimer situation</p>",collapse='')
  }
  
  htmlSource <- paste(htmlSource,"<p> <a href='#index'> <span class='backToIndexText'>Back to Index</span> </a> </p>",collapse='')
  }
  
  # -----------------------------
  # BARCODE RESULT
  # -----------------------------
  {
    htmlSource <- paste(htmlSource,"<a name='Barcodes'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h2> BARCODES </h2>",collapse='')
    
    # Show the definitions for assigned and filters and each of the columns
    htmlSource <- paste(htmlSource,"<a name='BarcodesDefinitions'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Definitions </h3>",collapse='')
    htmlSource <- paste(htmlSource,"<p> We have the following plots: </p>",collapse='')
    htmlSource <- paste(htmlSource,"<ul>
                                    <li> <b>Assigned:</b> We have four variables, how many assigned (grey) how many weren't reads (yellow) and how many where reads (green). </li>
                                    <li> <b>Filter:</b>  In here we summarize how many sequences we lost to the filter options, if they have N nucleotides (grey) if they didn't pass the minimum quality (red) if they didn't pass the minimum average quality (green) </li>
                                    </ul>",collapse='')
    htmlSource <- paste(htmlSource,"<p> We have the following columns: </p>",collapse='')
    htmlSource <- paste(htmlSource,"<ul>
                                    <li> <b>Barcode:</b> At vero eos et accusamus et iusto odio dignissimos </li>
                                    <li> <b>Total_Pre_N_Filter:</b> ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti </li>
                                    <li> <b>Total_Post_N_Filter:</b> quos dolores et quas molestias excepturi </li>
                                    <li> <b>Failed_N_Filter:</b> sint occaecati cupiditate non provident </li>
                                    <li> <b>Experiment_Unique_Sequences:</b> similique sunt in culpa qui officia deserunt mollitia animi </li>
                                    <li> <b>Barcode_Sum_Is_Sequence:</b> id est laborum et dolorum fuga. Et harum quidem rerum facilis est et expedita distinctio </li>
                                    <li> <b>Barcode_Sum_Total_Sequence:</b> Nam libero tempore, cum soluta nobis est eligendi optio cumque </li>
                                    <li> <b>Barcode_Sum_Is_Read:</b> nihil impedit quo minus id quod maxime placeat facere possimus </li>
                                    <li> <b>Barcode_Sum_Total_Reads:</b> omnis voluptas assumenda est, omnis dolor repellendus </li>
                                    <li> <b>Barcode_Unnasigned:</b> Temporibus autem quibusdam et aut officiis debitis aut rerum necessitatibus saepe eveniet </li>
                                    <li> <b>Barcode_non_read:</b> ut et voluptates repudiandae sint et molestiae non recusandae. Itaque earum rerum hic tenetur a sapiente delectus, ut aut reiciendis voluptatibus maiores alias consequatur aut perferendis doloribus asperiores repellat. </li>
                                    </ul>",collapse='')
      
    # Show the barcode_assigned plot.
    htmlSource <- paste(htmlSource,"<a name='BarcodesAssignedPlot'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Assigned Plot </h3>",collapse='')
    
    assignedPlotPath <- paste(plotFolderPath, "barcode_assigned.png" , sep="/")
    htmlSource <- paste(htmlSource,"<img src=",assignedPlotPath," alt='Assigned sequences for each barcode' height='20%' width='20%'>",collapse='')
    
    
    # Show the barcode_filter plot.
    htmlSource <- paste(htmlSource,"<a name='BarcodesFilterPlot'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Filter Plot </h3>",collapse='')
    
    filterPlotPath <- paste(plotFolderPath, "barcode_filters.png" , sep="/")
    htmlSource <- paste(htmlSource,"<img src=",filterPlotPath," alt='How many sequences passed the N filter for each barcode' height='20%' width='20%'>",collapse='')
  
    # Show an index with the barcode table result
    htmlSource <- paste(htmlSource,"<a name='BarcodesTable'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Barcode Table </h3>",collapse='')
    
    
    htmlSource <- paste(htmlSource , "<table>",collapse='')
    
    # --- Make the headers with the columns names
    barcodeColumnNames <- colnames(barcodeResultTable)
    htmlSource <- paste(htmlSource , "<tr>",collapse='')
    for(j in 1:ncol(barcodeResultTable)){
      htmlSource <- paste(htmlSource , "<th>",barcodeColumnNames[j],"</th>",collapse='')  
    }
    htmlSource <- paste(htmlSource , "</tr>",collapse='')
    
    # --- Make each of the rows
    for(i in 1:nrow(barcodeResultTable)){
      htmlSource <- paste(htmlSource , "<tr>",collapse='')  
      for(j in 1:ncol(barcodeResultTable)){
        
        htmlSource <- paste(htmlSource, "<td>",toString(barcodeResultTable[i,j]),"</td>", collapse='')
        
      }
      htmlSource <- paste(htmlSource , "</tr>",collapse='')
    }
    
    htmlSource <- paste(htmlSource , "</table>",collapse='')
    
    # The table has a function on_hover that tells you the definition of that column
  
    # For each of the index, show the individual plots of assigned and filter
    htmlSource <- paste(htmlSource,"<p> <a href='#index'> <span class='backToIndexText'>Back to Index</span> </a> </p>",collapse='')
  }
  # -----------------------------
  # CONFIG TABLE RESULT
  # -----------------------------
  {
    htmlSource <- paste(htmlSource,"<a name='Config'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h2> CONFIG TABLE </h2>",collapse='')
    htmlSource <- paste(htmlSource,"<a name='ConfigDefinitions'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Config definitions </h3>",collapse='')
    htmlSource <- paste(htmlSource,"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.",collapse='')
    # Show the definitions for cut rates
    # Show the definitions for frameshift
    # Show the definitions for the columns
    
    # Show the cut_rates plot.    
    htmlSource <- paste(htmlSource,"<a name='ConfigCutRatePlot'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Cut Rate Plot </h3>",collapse='')
    
    cutratePlotPath <- paste(plotFolderPath, "Cut_Rates.png" , sep="/")
    htmlSource <- paste(htmlSource,"<img src=",cutratePlotPath," alt='This is the Read Rate plot, also known as Cut Rate for inconsistances ' height='20%' width='20%'>",collapse='')
    
    
    # Show the frameshift plot.
    htmlSource <- paste(htmlSource,"<a name='ConfigFrameshiftPlot'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Frameshift Plot </h3>",collapse='')
    
    frameshiftPlotPath <- paste(plotFolderPath, "Frameshift.png" , sep="/")
    htmlSource <- paste(htmlSource,"<img src=",frameshiftPlotPath," alt='I'm trap in the server room, please send help!!' height='20%' width='20%'>",collapse='')
    
    # Show an index with the config table result
    htmlSource <- paste(htmlSource,"<a name='ConfigTable'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h3> Config Table </h3>",collapse='')
    
    htmlSource <- paste(htmlSource , "<table>",collapse='')
    
    # --- Make the headers with the columns names
    configColumnNames <- colnames(configResultTable)
    firstColumnName <- configColumnNames[1] # Save the name of the first column
    configColumnNames <- configColumnNames[2:length(configColumnNames)] # Delete the first column again
    configColumnNames <- c(firstColumnName, "archplot","histogram",configColumnNames) # Add two new columns at the beggining and append the first column once more

    # At the end we have Gene ID , archplot, histogram, the rest...
    

    
    
    htmlSource <- paste(htmlSource , "<tr>",collapse='')
    for(j in 1:length(configColumnNames)){
      htmlSource <- paste(htmlSource , "<th>",configColumnNames[j],"</th>",collapse='')  
    }
    htmlSource <- paste(htmlSource , "</tr>",collapse='')
    
    # --- Make each of the rows
    for(i in 1:nrow(configResultTable)){
      htmlSource <- paste(htmlSource , "<tr>",collapse='')  
      
      # First we add a name with a link to the alignments and a link to the unique table
      
      rowGeneID <- toString(configResultTable[i,1])
      rowBarcode <- toString(configResultTable[i,2])
      
      htmlSource <- paste(htmlSource, '<td> <div class= "IDContainer"> <div class="geneID"> <p>',rowGeneID,'</p> </div> <div class="alignmentLink"> <a href="',htmlFolder,"/",rowGeneID, '_', rowBarcode,'_alignments.html"> <span class="configTableAlignmentLink"> alignments </span> </a> </div> <div class="uniquesLink"> <a href="#',rowGeneID,'_unique"> <span class="configTableAlignmentLink"> details </span> </a> </div> </div> </td>', collapse='', sep='')
            
      # Then we add two columns, one with the arch plot and another with the histogram
      # For each of the index, show the individual plots of abundance and archplot
      
      currentPlot <- paste(plotFolderPath, "/", rowGeneID, "_" , rowBarcode , "_archplot.png" , sep='', collapse='')
      if(file.exists(currentPlot)){
        
        #Set the name for the new PNG file
        targetPNGFile <- paste(htmlFolder,"/",rowGeneID, '_', rowBarcode,"_archplot.png", sep='', collapse='')
        print(targetPNGFile)
        
        #Read the original image
        archPlotRaster <- readPNG(currentPlot)
        
        
        # Save the same file but 10 times smaller with some AA
        png(targetPNGFile, width=1181, height=210, units='px')
        par(mar = c(0,0,0,0) ) # No margins around the image
        plot(c(0, 1181), c(0, 210), type = "n", xlab = "", ylab = "")
        Axis(side=1, labels=FALSE)
        Axis(side=2, labels=FALSE)
        rasterImage(archPlotRaster, 0, 0, 1181, 210, interpolate = TRUE)
        dev.off()
        
        # Show the graphic in the table with a link to the original       
        htmlSource <- paste(htmlSource , "<td class='archPlotCell'> <a href=",currentPlot,"> <img src=",targetPNGFile," alt='Archplot, blue rectangle is the AP window' >  </a> </td>",collapse='')  
      }
      else{
        htmlSource <- paste(htmlSource , "<td> Archplot not found  </td>",collapse='')  
      }
            
      currentPlot <- paste(plotFolderPath, "/", rowGeneID, "_" , rowBarcode , ".png" , sep='', collapse='')
      if(file.exists(currentPlot)){
        htmlSource <- paste(htmlSource , "<td> <img src=",currentPlot," alt='Showing 15 or less (all)' height='10%'> </td>",collapse='')  
      }
      else{
        htmlSource <- paste(htmlSource , "<td> Histogram not found  </td>",collapse='')  
      }
      
      
      for(j in 2:ncol(configResultTable)){ # Add the rest of the row
        
        htmlSource <- paste(htmlSource, "<td>",toString(configResultTable[i,j]),"</td>", collapse='')
        
      }
      htmlSource <- paste(htmlSource , "</tr>",collapse='')
    }
    
    htmlSource <- paste(htmlSource , "</table>",collapse='')
    
    
    
    # The table has a function on_hover that tells you the definition of that column
    
    # The index will link inside the page to the place where the unique table for each ID is located
    
    
    htmlSource <- paste(htmlSource,"<p> <a href='#index'> <span class='backToIndexText'>Back to Index</span> </a> </p>",collapse='')
  
  }
  # -----------------------------
  # LOG RESULT
  # -----------------------------
  {
    htmlSource <- paste(htmlSource,"<a name='Log'></a> ",collapse='') 
    htmlSource <- paste(htmlSource,"<h2> LOGS </h2>",collapse='')
    htmlSource <- paste(htmlSource,"Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi consequatur? Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam nihil molestiae consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?",collapse='')
    
    htmlSource <- paste(htmlSource,"<p> <a href='#index'> <span class='backToIndexText'>Back to Index</span> </a> </p>",collapse='')
    
    
  }
  # -----------------------------
  # ALIGNMENTS RESULT
  # -----------------------------
  {
    htmlSource <- paste(htmlSource,"<a name='Barcodes'></a> ",collapse='')
    htmlSource <- paste(htmlSource,"<h2> Alignments </h2>",collapse='')  
    
    # Make each of the alignments summary
    for(i in 1:nrow(configResultTable)){
      
      #Locate the file
      rowGeneID <- toString(configResultTable[i,1])
      rowBarcode <- toString(configResultTable[i,2])
      currentAligmentPath <- paste(resultFolderPath, "/", rowGeneID, "_" , rowBarcode , "/alignments.txt" , sep='', collapse='')
      
      # Make a new HTML file with this information
      htmlAlignmentFilePath <- paste(htmlFolder,"/",rowGeneID, "_" , rowBarcode, "_alignments.html",sep='')            
      htmlAlignmentFD<-file(htmlAlignmentFilePath)
      
      #print(currentAligmentPath)

      htmlAlignment <-" <!DOCTYPE html> <html> <head> <style> </style> </head> <body>"
      
      htmlAlignment <- paste(htmlAlignment,'<a name="',rowGeneID,'_', rowBarcode,'_alignments"></a>',collapse='' , sep='')
      
      # If the master alignment for that combo exist show it, otherwise tell that you couldn't find it
      if(file.exists(currentAligmentPath)){
      
        alignmentFD <-  file(currentAligmentPath, open="r")
        alignmentText <- readLines(currentAligmentPath,encoding="UTF-8")
        close(alignmentFD)
        
        htmlAlignment <- paste(htmlAlignment,"<p>",alignmentText,"</p>",collapse='')
        #htmlAlignment <- paste(htmlAlignment,"<p> File found</p>",collapse='')
        
      }
      else{
        
        htmlAlignment <- paste(htmlAlignment,"<p> Alignment summary not found</p>",collapse='')
      
      }
      
      htmlAlignment <- paste(htmlAlignment,"</body> </html>")
      
      # Write into file and close the file descriptor
      writeLines(htmlAlignment, htmlAlignmentFD)
      close(htmlAlignmentFD)
      
      
    }
    
  }



  # Close the headers and any other <> that need clousure
  htmlSource <- paste(htmlSource, htmlEnding, collapse='')
  
  # Write the HTML code in disk; inside the HTML folder


  htmlFilePath <- paste(htmlFolder,"/report.html",sep='')            
  htmlSourceFD<-file(htmlFilePath)
  writeLines(htmlSource, htmlSourceFD)
  close(htmlSourceFD)
  
  

}