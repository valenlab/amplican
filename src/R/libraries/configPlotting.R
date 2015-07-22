makePlots <- function(CUT_RATE, NO_CUT_RATE, NO_DELETION, IN_FRAME, FRAMESHIFT, MULTIPLES,
                      DELETION_COLOR, CUT_COLOR, POSITION_WINDOW, START, END, DUPLICATE_INTO_FOLDERS,
                      PLOTBACK, PLOTGRID, LINETYPE, FREQONE, ACCUMULATE, FILTERENDS, SIZE,
                      analysisDataframe, processID, ANALYSIS_PATH, plotPath, dataPath, alignmentDistance, FORMAT){
  
  ALIGNMENT_DISTANCE <- alignmentDistance
  
  # Get some timing variables ready
  {
    t1 <- Sys.time()
    t2 <- Sys.time()
    t3 <- Sys.time()
    t4 <- Sys.time()
    td <- as.numeric(t2-t1, units = "secs")
    td2 <- as.numeric(t4-t3, units = "secs")
    
    histogramPlottingTiming <- 0
    archplotTiming          <- 0
  }
  
  print(paste("Nice to meet you; I'm mighty processor number: ",processID))
  
  # For each line of the config table
  for (j in 1:nrow(analysisDataframe)){
  #for (j in 1:1){
  
    print(paste("The processor",processID,"is doing line",j,"of",nrow(analysisDataframe), "Progress of: ",round(j/nrow(analysisDataframe)*100,2),"%"))
    
    #------------------------------------------------------
    # GET BASIC INFO
    #------------------------------------------------------
    {

      # Set the time to zero
      histogramPlottingTiming <- 0
      archplotTiming          <- 0
      
      # Get the name of the row, that will tell us the folder
      currentID <- analysisDataframe[j,"ID"]
      currentBarcode <- analysisDataframe[j,"Barcode"]
      
      # Get the amplicon
      amplicon <- analysisDataframe[j,"Genome"]
      
      # Get the position of the alignment and how wide it is
      # TODO: Adjust for more than one
      allPositions <- countUppercaseGroups(toString(amplicon))
      alignmentPositions <- allPositions[[2]][1]
      widePosition <- allPositions[[3]][1]
      
      # From that folder read the unique table and the cut table
      dataFolder <- paste(ANALYSIS_PATH,"/",currentID,"_",currentBarcode,sep='')
      
      uniquesDF  <- read.table(paste(dataFolder, "/", currentID,"_",currentBarcode,"_uniques.txt",   sep=''), stringsAsFactors=FALSE, header=TRUE, sep="\t")
      cutsDF     <- read.table(paste(dataFolder, "/", currentID,"_",currentBarcode,"_cutsTable.txt", sep=''), stringsAsFactors=FALSE, header=TRUE, sep="\t")
      
    }

    #------------------------------------------------------
    # PLOT THE HISTOGRAMS
    #------------------------------------------------------
    {
      t1 <- Sys.time()
      
      uniqueGraphics <- data.frame()
      
      # TOP 15 GOOD READS
      #-------------------
      
      # --  If we don't have an uniqueDF empty, make a subset with the reads and sorted by total
      if(nrow(uniquesDF) >= 1){
        uniqueGraphics <- uniquesDF[,c("Total","Is_Read", "Alignment", "Forward", "Reverse")]        
        uniqueGraphics <- subset(uniqueGraphics, Is_Read == TRUE)
        uniqueGraphics <- uniqueGraphics[order(-uniqueGraphics$Total), ]
      }
      
      # -- Now we need to sort the x axis, we are going to represent the name of the alignment there
      # -- but we don't want the name of the alignment as the sort. So we are going to use the "rowname" plus
      # -- the name of the alignment to trick ggplot into sort it like that
      
      # TODO: Fix this flow
      if(nrow(uniquesDF) >= 1){
        
        if(nrow(uniqueGraphics) >= 1){
        
          baseZero <- nchar(toString(nrow(uniqueGraphics)))
          
          for(p in 1:nrow(uniqueGraphics)){
            
            uniqueGraphics$Alignment[p] <- rownames(uniqueGraphics[p,])
            
            #Find how many zeros we need to fix the sort
            totalZeros  <- baseZero - nchar(toString(p))
            stringZeros <- paste(rep("0",totalZeros),collapse='')
            
            uniqueGraphics$Alignment[p] <- paste(stringZeros,p,"_*_",toString(uniqueGraphics$Alignment[p]),sep='')
            
          }
          
          selectedGraphics <- NULL
          
          # Get maximum 15 rows, we only care about the first rows that have the highest total
          if(nrow(uniqueGraphics) > 15){
            selectedGraphics <- uniqueGraphics[1:15,]##ggplot is complaining when I do this on the fly on the next line
          }
          else{
            selectedGraphics <- uniqueGraphics[1:nrow(uniqueGraphics),]
          }
          
          # Change the names of the columns
          # Alignments -> Sequence
          # Total      -> Reads
          colnames(selectedGraphics) <- c("Reads", "Is_Cut", "Sequence", "Forward", "Reverse")
          
          histogramPlot <- ggplot(data=selectedGraphics, aes(x=Sequence, y=Reads), na.rm = TRUE)
          histogramPlot <- histogramPlot + geom_bar(stat="identity", na.rm = TRUE, fill = CUT_RATE)          
          histogramPlot <- histogramPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                 panel.background = element_rect(fill = PLOTBACK),
                                                 panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
                                                 panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
                                                 panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
          
          # Save the file into the plot folder
          ggsave(paste(plotPath,"/",currentID,"_",currentBarcode,"_histogramGOODCUT.", FORMAT, sep = ''))
          
          fileName <- paste(dataPath,"/",currentID,"_",currentBarcode,"_histogramGOODCUT.txt",sep = '', collapse = '')
          write.table(selectedGraphics, file = fileName, quote = FALSE, sep = "\t")
          
        }  
        
      }
      
      
      # TOP 15 BAD READS
      # ----------------------
      
      # --  If we don't have an uniqueDF empty, make a subset with the reads and sorted by total
      if(nrow(uniquesDF) >= 1){
        uniqueGraphics <- uniquesDF[,c("Total","Is_Read", "Alignment", "Forward", "Reverse")]        
        uniqueGraphics <- subset(uniqueGraphics, Is_Read == FALSE)
        uniqueGraphics <- uniqueGraphics[order(-uniqueGraphics$Total), ]
      }
      
      # -- Now we need to sort the x axis, we are going to represent the name of the alignment there
      # -- but we don't want the name of the alignment as the sort. So we are going to use the "rowname" plus
      # -- the name of the alignment to trick ggplot into sort it like that
      
      # TODO: Fix this flow
      if(nrow(uniquesDF) >= 1){
        
        if(nrow(uniqueGraphics) >= 1){
          
          baseZero <- nchar(toString(nrow(uniqueGraphics)))
          
          for(p in 1:nrow(uniqueGraphics)){
            
            uniqueGraphics$Alignment[p] <- rownames(uniqueGraphics[p,])
            
            #Find how many zeros we need to fix the sort
            totalZeros  <- baseZero - nchar(toString(p))
            stringZeros <- paste(rep("0",totalZeros),collapse='')
            
            uniqueGraphics$Alignment[p] <- paste(stringZeros,p,"_*_",toString(uniqueGraphics$Alignment[p]),sep='')
            
          }
          
          selectedGraphics <- NULL
          
          # Get maximum 15 rows, we only care about the first rows that have the highest total
          if(nrow(uniqueGraphics) > 15){
            selectedGraphics <- uniqueGraphics[1:15,]##ggplot is complaining when I do this on the fly on the next line
          }
          else{
            selectedGraphics <- uniqueGraphics[1:nrow(uniqueGraphics),]
          }
          
          # Change the names of the columns
          # Alignments -> Sequence
          # Total      -> Reads
          colnames(selectedGraphics) <- c("Reads", "Is_Cut", "Sequence", "Forward", "Reverse")
          
          #         ggplot(data=selectedGraphics, aes(x=Alignment, y=Total), na.rm = TRUE)
          #               + geom_bar(stat="identity", na.rm = TRUE)
          #               + theme(axis.text.x = element_text(angle = 90, hjust = 1))
          
          histogramPlot <- ggplot(data=selectedGraphics, aes(x=Sequence, y=Reads), na.rm = TRUE)
          histogramPlot <- histogramPlot + geom_bar(stat="identity", na.rm = TRUE, fill = NO_CUT_RATE)
          histogramPlot <- histogramPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                 panel.background = element_rect(fill = PLOTBACK),
                                                 panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
                                                 panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
                                                 panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
          
          # Save the file into the plot folder
          ggsave(paste(plotPath,"/",currentID,"_",currentBarcode,"_histogramNOCUT.", FORMAT, sep = ''))
          
          fileName <- paste(dataPath,"/",currentID,"_",currentBarcode,"_histogramNOCUT.txt",sep = '', collapse = '')
          write.table(selectedGraphics, file = fileName, quote = FALSE, sep = "\t")
          
        }  
        
      }
      
      
      # TOP 15 GOOD OR BAD READS
      # ---------------------------
      
      # --  If we don't have an uniqueDF empty, make a subset with the reads and sorted by total
      if(nrow(uniquesDF) >= 1){
        uniqueGraphics <- uniquesDF[,c("Total","Is_Read", "Alignment", "Forward", "Reverse")]        
        uniqueGraphics <- uniqueGraphics[order(-uniqueGraphics$Total), ]
      }
      
      # -- Now we need to sort the x axis, we are going to represent the name of the alignment there
      # -- but we don't want the name of the alignment as the sort. So we are going to use the "rowname" plus
      # -- the name of the alignment to trick ggplot into sort it like that
      
      # TODO: Fix this flow
      if(nrow(uniquesDF) >= 1){
        
        if(nrow(uniqueGraphics) >= 1){
          
          baseZero <- nchar(toString(nrow(uniqueGraphics)))
          
          for(p in 1:nrow(uniqueGraphics)){
            
            uniqueGraphics$Alignment[p] <- rownames(uniqueGraphics[p,])
            
            #Find how many zeros we need to fix the sort
            totalZeros  <- baseZero - nchar(toString(p))
            stringZeros <- paste(rep("0",totalZeros),collapse='')
            
            uniqueGraphics$Alignment[p] <- paste(stringZeros,p,"_*_",toString(uniqueGraphics$Alignment[p]),sep='')
            
          }
          
          selectedGraphics <- NULL
          
          # Get maximum 15 rows, we only care about the first rows that have the highest total
          if(nrow(uniqueGraphics) > 15){
            selectedGraphics <- uniqueGraphics[1:15,]##ggplot is complaining when I do this on the fly on the next line
          }
          else{
            selectedGraphics <- uniqueGraphics[1:nrow(uniqueGraphics),]
          }
          
          # Change the names of the columns
          # Alignments -> Sequence
          # Total      -> Reads
          colnames(selectedGraphics) <- c("Reads", "Is_Cut", "Sequence", "Forward", "Reverse")
          
          #         ggplot(data=selectedGraphics, aes(x=Alignment, y=Total), na.rm = TRUE)
          #               + geom_bar(stat="identity", na.rm = TRUE)
          #               + theme(axis.text.x = element_text(angle = 90, hjust = 1))
          
          histogramPlot <- ggplot(data=selectedGraphics, aes(x=Sequence, y=Reads, fill = Is_Cut), na.rm = TRUE)
          histogramPlot <- histogramPlot + geom_bar(stat="identity", na.rm = TRUE)
          histogramPlot <- histogramPlot + scale_colour_discrete(drop=TRUE , limits = c("FALSE","TRUE"))
          histogramPlot <- histogramPlot + scale_fill_manual(values=c(NO_CUT_RATE, CUT_RATE))
          histogramPlot <- histogramPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                 panel.background = element_rect(fill = PLOTBACK),
                                                 panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
                                                 panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
                                                 panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
          
          # Save the file into the plot folder
          ggsave(paste(plotPath,"/",currentID,"_",currentBarcode,"_histogramBOTHCUT.", FORMAT, sep = ''))
          
          fileName <- paste(dataPath,"/",currentID,"_",currentBarcode,"_histogramBOTHCUT.txt",sep = '', collapse = '')
          write.table(selectedGraphics, file = fileName, quote = FALSE, sep = "\t")
          
        }  
        
      }
      
      
      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      histogramPlottingTiming <- td
  
    }
    
    #------------------------------------------------------
    # PLOT THE RESULTS
    # -- Plot the archplot
    # -- Plot the mutations
    #------------------------------------------------------
    {
      
      t1 <- Sys.time()
      
      # Here are the values for the different configurations of SIZE. The first column is the font size
      # the second column is the width for the resulting image. The last one is the SIZE value (string)
      # 1425 ; 1000 = 4K
      #      = 1080
      #      ; 62   = 720
      #      = VGA
      
      baseWidth   <- 0
      baseSize    <- 0
      legendSize  <- 0
      legendTitle <- 0
      
      if(SIZE == "4K"){
        baseSize    <- 1425
        baseWidth   <- 1000
        legendSize  <- 6
        legendTitle <- 10
      }
      else{
        
        if(SIZE == "1080"){
          baseSize    <- 1425
          baseWidth   <- 1000
          legendSize  <- 6
          legendTitle <- 10
        }
        else{
          
          if(SIZE == "720"){
            baseSize    <- 89
            baseWidth   <- 62
            legendSize  <- 2
            legendTitle <- 3
          }
          else{
            baseSize    <- 1425
            baseWidth   <- 1000
            legendSize  <- 6
            legendTitle <- 10
          }
        }
      }
      
      # Lets find out the length of the amplicon first and generate a special DF
      # This DF is use to plot the nucleotides inside the archplot and the mutation plots
      # DONT try to do this with a for because it doubles the time (trust me).
      
      ampliconLength <- nchar(toString(amplicon)) # We need this variable here
      addjustedSize  <- baseSize/ampliconLength  # TODO: Not sure that this is correct yet
      
      ampliconDF <- data.frame(matrix(NA,ampliconLength,2))
      colnames(ampliconDF) <- c("Letter", "X")
      ampliconDF$Letter <- strsplit(toString(amplicon),"")[[1]]
      ampliconDF$X      <- seq(1,ampliconLength)
      
      
   
      
      # Now lets make the archplot
      
      # First lets delete the bad rows that comes from the trace of the forward or reverse sequence.
      # This step is totally optional and up to the user.
      if(FILTERENDS == TRUE && nrow(cutsDF)>0){
        cutsDF <- subset(cutsDF, !(Type=="deletion" & Start == 1) & !(Type=="deletion" & End == ampliconLength))
      }
      
      if(nrow(cutsDF)>0){
        
        # First lets construct the GRanges data
        minStart       <- min(0,min(cutsDF$Start))
        maxEnd         <- max(nchar(toString(amplicon)), max(cutsDF$End))
        graphicWide    <- maxEnd - minStart
        maxFreq        <- 1
        startOffset    <- minStart*(-1)
        endOffset      <- 0
        
        if(FREQONE==FALSE){
          maxFreq      <- max(cutsDF$Freq)
        }
        
        if(maxEnd>ampliconLength){
          endOffset <- max(cutsDF$End) - ampliconLength
        }
        
        addjustedSize  <- 1425/(ampliconLength + sqrt(startOffset*startOffset) + endOffset) #1425 is manually set up for a image size of 11K pixels. Should be a linear proportion there
        
        # We need to accumulate the cuts that start and end at the same position, IF the user want to do it.
        if(ACCUMULATE == TRUE){
          cdata         <- ddply(cutsDF, c("Start", "End", "Type", "Location"), summarise, Freq = sum(Freq))
          cdata$Wide    <- cdata$End - cdata$Start
          cdata$Gene_ID <- cutsDF$Gene_ID[1]
          cutsDF        <- cdata
        }        
        
        
        gr <- GRanges(
          seqnames = cutsDF$Gene_ID,
          
          IRanges(
            start  = cutsDF$Start,
            width  = cutsDF$Wide,
          ),
          value  = cutsDF$Freq,
          Type   = cutsDF$Type,
          Frequency_Relative = sqrt(cutsDF$Freq*cutsDF$Freq)

        )

        # And now, plot it and save it to disk
        myArchplot <- ggplot(gr) + #Plot the Granges data
          #In an archplot, with height and size proportional to the frequency
          geom_arch( aes(color = Type, height = value, size = Frequency_Relative ), alpha = 0.5) + 
          # Keep leyend consistant between plot. Red is a deletion, blue is a cut always.
          scale_colour_discrete(drop=TRUE , limits = c("deletion","cut")) + 
          #scale_colour_manual(values = c(CUT_RATE, NO_CUT_RATE)) + # This will trigger an adding another scale, but it seems to be the only way! NOTWORKING NOW!!!??? WTF R!!!
          scale_fill_manual(values=c(NO_CUT_RATE, CUT_RATE)) +
          # Make a green line at the end of the amplicon, +1 makes the line after the last letter of the amplicon
          annotate("pointrange", x = nchar(toString(amplicon)) + 1 , y = 0, ymin = -maxFreq, ymax = maxFreq, colour = END, size = 1) + 
          # Make a green line at the start of the amplicon       
          annotate("pointrange", x = 1, y = 0, ymin = -maxFreq, ymax = maxFreq, colour = START, size = 1) +
          # Make the bw background
          theme(
            panel.background = element_rect(fill = PLOTBACK),
            panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
            panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
            panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
        

        # For each of the alignment position (usually only one), make a rectangle with the AP position, and the range around it.
        if(alignmentPositions>=1){
          for(k in 1:alignmentPositions){
            
            myArchplot <- myArchplot + annotate("rect", xmin = alignmentPositions[k], xmax = alignmentPositions[k]+widePosition[k], ymin = -maxFreq/5, ymax = maxFreq/5, alpha = .2) 
            myArchplot <- myArchplot + annotate("rect", xmin = alignmentPositions[k]-ALIGNMENT_DISTANCE, xmax = alignmentPositions[k]+widePosition[k]+ALIGNMENT_DISTANCE, ymin = -maxFreq/5, ymax = maxFreq/5, alpha = .1, colour = POSITION_WINDOW) 
            
          }}

        myArchplot <- myArchplot + geom_text(data = ampliconDF, x = ampliconDF$X, y = 0, hjust = 0, label = ampliconDF$Letter, size=addjustedSize, aes(family="mono")) 
        
        # For each of the characters on the amplicon, lets write one letter
#         for(k in 1:ampliconLength){
# 
#           myArchplot <- myArchplot + geom_text(data = NULL, x = k, y = 0, hjust = 0, label = strsplit(toString(amplicon),"")[[1]][k] , size=addjustedSize, aes(family="mono")) 
#           
#         }

        # Save the file into the plot folder
        ggsave(paste(plotPath,"/",currentID,"_",currentBarcode,"_archplot.", FORMAT ,sep = ''), width = 1000 , units="mm", limitsize=FALSE)
        
        fileName <- paste(dataPath,"/",currentID,"_",currentBarcode,"_archplot.txt",sep = '', collapse = '')
        write.table(cutsDF, file = fileName, quote = FALSE, sep = "\t")
        
      }
      
      
      # Now lets make the mutation plot.
      # We need the mutation data from the unique table for both forward and reverse
      # In order to do that we are goint to fill this two dataframes of size =  4 x amplicon length.
      
      # The first  dataframe counts the mutations in the forward sequences
      # The second dataframe counts the mutations in the reverse sequences
      # Both, only counts mutations from the amplicon point of view.
      
      # We have 4 possible mutations. From X to A,T,C,G (one of them will always be 0)
      
      forwardMutations <- data.frame(matrix(0,ampliconLength,4))
      reverseMutations <- data.frame(matrix(0,ampliconLength,4))
      colnames(forwardMutations) <- c("A", "T", "C", "G")
      colnames(reverseMutations) <- c("A", "T", "C", "G")
      
      # Now we need to go through the whole unique matrix and count the mutations for each
      
      # Useful variables: "Total"  "Forward_Deletions_Genome_Coordinates_String"	"Reverse_Deletions_Genome_Coordinates_String"			"Is_Read"	
      if(nrow(uniquesDF)>0){
      for(i in 1:nrow(uniquesDF)){
        
        # Get the info from the unique table
        forwardAmpliconRelativeLiteString <- as.character(uniquesDF$Forward_Deletions_Genome_Coordinates_String[i])
        reverseAmpliconRelativeLiteString <- as.character(uniquesDF$Reverse_Deletions_Genome_Coordinates_String[i])
        
        #@0@!@2@47,53*158,202*!@1@77,t,C*"
        #alignmentMutationsOriginal, alignmentMutationsMutated, alignmentMutationsPosition
        
        # Transform the info in data
        forwardAmpliconData             <- getEventInfo(forwardAmpliconRelativeLiteString)
        reverseAmpliconData             <- getEventInfo(reverseAmpliconRelativeLiteString)
        
        #forwardGenomeMutationsOriginals   <- forwardGenomeData[[7]]
        forwardAmpliconMutationsMutated     <- forwardAmpliconData[[8]]
        forwardAmpliconMutationsPositions   <- forwardAmpliconData[[9]]
        
        #reverseGenomeMutationsOriginals   <- reverseGenomeData[[7]]
        reverseAmpliconMutationsMutated     <- reverseAmpliconData[[8]]
        reverseAmpliconMutationsPositions   <- reverseAmpliconData[[9]]
        
        # Add the data to our forward result matrix
        if(length(forwardAmpliconMutationsMutated)>0){
        for(k in 1:length(forwardAmpliconMutationsMutated)){
          
          # Initialize the proper variables
          column <- 1
          forwardAmpliconMutationsMutated[k] <- toupper(forwardAmpliconMutationsMutated[k])
          
          # Select the appropiate column
          if(     forwardAmpliconMutationsMutated[k] == "T"){ column <- 2}
          else{if(forwardAmpliconMutationsMutated[k] == "C"){ column <- 3}
          else{if(forwardAmpliconMutationsMutated[k] == "G"){ column <- 4}}}
        
          # Add the result to our dataframes
          forwardMutations[forwardAmpliconMutationsPositions[k],column] <- forwardMutations[forwardAmpliconMutationsPositions[k],column] + uniquesDF$Total[i]
        
        }}
        
        # Add the data to our reverse result matrix
        if(length(reverseAmpliconMutationsMutated)>0){
          for(k in 1:length(reverseAmpliconMutationsMutated)){
            
            # Initialize the proper variables
            column <- 1
            reverseAmpliconMutationsMutated[k] <- toupper(reverseAmpliconMutationsMutated[k])
            
            # Select the appropiate column
            if(     reverseAmpliconMutationsMutated[k] == "T"){ column <- 2}
            else{if(reverseAmpliconMutationsMutated[k] == "C"){ column <- 3}
            else{if(reverseAmpliconMutationsMutated[k] == "G"){ column <- 4}}}
            
            # Add the result to our dataframes
            reverseMutations[reverseAmpliconMutationsPositions[k],column] <- reverseMutations[reverseAmpliconMutationsPositions[k],column] - uniquesDF$Total[i]
            
         }}
      
      
      }}

      # Write down the data for the mutation without the melting, so just as it is right now
      fileName <- paste(dataPath,"/",currentID,"_",currentBarcode,"_forwardMutation.txt",sep = '', collapse = '')
      write.table(forwardMutations, file = fileName, quote = FALSE, sep = "\t")

      fileName <- paste(dataPath,"/",currentID,"_",currentBarcode,"_reverseMutation.txt",sep = '', collapse = '')
      write.table(reverseMutations, file = fileName, quote = FALSE, sep = "\t")

      # Add the position column and melt the data
      forwardMutations$Position <- row.names(forwardMutations)
      reverseMutations$Position <- row.names(reverseMutations)

      # Melt everything
      forwardMelted <- melt(forwardMutations, id.var="Position")
      reverseMelted <- melt(reverseMutations, id.var="Position")

      # Change the columns names
      colnames(forwardMelted) <- c("Position", "Nucleotide","Frequency")
      colnames(reverseMelted) <- c("Position", "Nucleotide","Frequency")

      # Find out the maximum frequency, so the plot is centered at 0
      maxFreq      <- max(forwardMelted$Frequency)
      minFreq      <- min(reverseMelted$Frequency)
      if(maxFreq < minFreq * (-1)) {maxFreq <- minFreq * (-1)}

      # Now plot the mutations and save it to disk


#       mutationPlot <- ggplot(data = bothMelted, aes(x = as.numeric(Position) + 0.5, y = Frequency, fill = Nucleotide)) +
        mutationPlot <- ggplot() + 
                        geom_bar(data = forwardMelted, aes(x=as.numeric(Position) + 0.5, y=Frequency, fill=Nucleotide),stat = "identity", binwidth = 1) +
                        geom_bar(data = reverseMelted, aes(x=as.numeric(Position) + 0.5, y=Frequency, fill=Nucleotide),stat = "identity", binwidth = 1) +
#                       geom_bar(stat = "identity", binwidth = 1) + 
                      # Make a green line at the end of the amplicon, +1 makes the line after the last letter of the amplicon
                      annotate("pointrange", x = nchar(toString(amplicon)) + 1 , y = 0, ymin = -maxFreq, ymax = maxFreq, colour = END, size = 1) + 
                      # Make a green line at the start of the amplicon       
                      annotate("pointrange", x = 1, y = 0, ymin = -maxFreq, ymax = maxFreq, colour = START, size = 1) +
                      # Make the bw background
                      theme(
                        panel.background = element_rect(fill = PLOTBACK),
                        panel.grid.major = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.2),
                        panel.grid.minor = element_line(colour = PLOTGRID, linetype = LINETYPE, size = 0.5),
                        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

      # For each of the alignment position (usually only one), make a rectangle with the AP position, and the range around it.
      if(alignmentPositions>=1){
      for(k in 1:alignmentPositions){
        
        mutationPlot <- mutationPlot + annotate("rect", xmin = alignmentPositions[k], xmax = alignmentPositions[k]+widePosition[k], ymin = -maxFreq/5, ymax = maxFreq/5, alpha = .2) 
        mutationPlot <- mutationPlot + annotate("rect", xmin = alignmentPositions[k]-ALIGNMENT_DISTANCE, xmax = alignmentPositions[k]+widePosition[k]+ALIGNMENT_DISTANCE, ymin = -maxFreq/5, ymax = maxFreq/5, alpha = .1, colour = POSITION_WINDOW) 
                          
      }}

      mutationPlot <- mutationPlot + geom_text(data = ampliconDF, x = ampliconDF$X, y = 0, hjust = 0, label = ampliconDF$Letter, size=addjustedSize, aes(family="mono")) 

      # For each of the characters on the amplicon, lets write one letter
#       for(k in 1:ampliconLength){
#       
#         mutationPlot <- mutationPlot + geom_text(data = NULL, x = k, y = 0, hjust = 0, label = strsplit(toString(amplicon),"")[[1]][k] , size=addjustedSize, aes(family="mono")) 
#                         
#       }
  
      # Save the file into the plot folder
      ggsave(paste(plotPath,"/",currentID,"_",currentBarcode,"_mutations.", FORMAT ,sep = ''), width = 1000 , units="mm", limitsize=FALSE)      

      t2 <- Sys.time()
      td <- as.numeric(t2-t1, units = "secs")
      archplotTiming <- td
      
    }
    
    
  }
  
  print(paste(processID, "finish!"))

}


