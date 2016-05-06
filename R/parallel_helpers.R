#' This function takes a dataframe and a number of processors and gives back the
#' same dataframe divided into several subsets. Every subset crafted so each
#' processors takes aproximatly the same amount of rows. Note that that doesn't
#' mean that they will take the same workload. This is just a arbitrary way of
#' dividing the work.
#'
#' See:
#'
#'   divideWorkBySize()
#'   divideWorkByBarcode()
#'@param   totalProcessors: (Int)       An integer with the number of processors.
#'                                This is also the amount of groups that will
#'                                result from this function. The number must be
#'                                1 or bigger.
#'@param   configDataframe: (Dataframe) A dataframe you want to divide. It doesn't
#'                                require any special column. But it must not be
#'                                empty. It can be full of NA rows thought.
#'@return   (List) A list with several dataframes. Each dataframe is a subgroup of
#'          configDataframe. If the number of rows of configDataframe is bigger
#'          than totalProcessors, the list have length totalProcessors. If it is
#'          smaller, the length is equal to the number of rows of
#'          configDataframe.
divideWork <- function(totalProcessors, configDataframe){
  # Make some initial calculations
  totalConfigLines <- nrow(configDataframe)
  processorsLines <- rep(as.integer(totalConfigLines/totalProcessors),totalProcessors)

  if((totalConfigLines%%totalProcessors)>0){
    for(i in 1:((totalConfigLines%%totalProcessors)) ){
      processorsLines[i] <- processorsLines[i] + 1
    }
  }

  processorsSubDataframes <- list()

  start <- 1
  end <- 0
  offset <- 0

  for(i in 1:totalProcessors){
    # Find out the ending line
    end <- offset + processorsLines[i]
    offset <- end
    #     print(paste("For processor",i))
    #     print(start)
    #     print(end)
    # Make a subset based on that
    subData <- NULL
    if(start <= end){
      subData <- configDataframe[start:end,]
    }
    processorsSubDataframes[[i]] <- subData
    # Update the next start to the current end + 1
    start <- end + 1
  }
  return(processorsSubDataframes)
}
#' This function divide the dataframe into similar dataframes based on the
#' colName values. The function try to divide the work so each processor gets
#' aproximaly the same amount of colName summations. This is commonly known as
#' the PARTITION PROBLEM. This implementation uses a greedy approach with order
#' O(nlogn).
#'
#' See:
#'
#'   divideWork()
#'   divideWorkByBarcode()
#'
#'@param   totalProcessors: (Int)       An integer with the number of processors.
#'                                This is also the amount of groups that will
#'                                result from this function. The number must be
#'                                1 or bigger.
#'
#'@param   configDataframe: (Dataframe) A dataframe you want to divide. It must not be
#'                                empty.
#'
#'@param   colName:         (String)    The name of a column in configDataframe. This
#'                                column should be full of of values that can be
#'                                sorted and can be added (typically either int
#'                                or float). It can be negative values.
#'@return   (List) A list with several dataframes. Each dataframe is a subgroup of
#'          configDataframe. If the number of rows of configDataframe is bigger
#'          than totalProcessors, the list have length totalProcessors. If it is
#'          smaller, the length is equal to the number of rows of
#'          configDataframe.
divideWorkBySize <- function(totalProcessors, configDataframe, colName){
  # In here we are going to store how many lines goes for each processor
  # at the beggining is a list of empty lists (slow)
  processorsLines <- rep( list(list()), totalProcessors )

  # In here we store the sum for each processor
  processorsSums <- rep(0, totalProcessors)

  # Get the data sorted by value in decreased order
  rowsSorted <- order(configDataframe[,c(colName)], decreasing = TRUE)

  # Get the total of lines
  totalConfigLines <- nrow(configDataframe)

  # Initialize the lists. Each one get the top rows
  for(i in 1:totalProcessors){
    processorsLines[i] <- c(rowsSorted[i])
    processorsSums[i]  <- configDataframe[,c(colName)][rowsSorted[i]]
  }

  # Now; for the rest of the non-top rows, added to set that will render the smaller sum
  if(totalConfigLines > 1){
    for(i in (totalProcessors+1):totalConfigLines){

      #Initialize the variables for this iteration
      analyzeLine  <- rowsSorted[i]
      currentValue <- configDataframe[,c(colName)][analyzeLine]

      # Minimum and reference
      minimum <- processorsSums[1] + currentValue
      reference <- 1

      # Try all combinations (except the first one which is minimum by default)
      if(totalProcessors>=2){
        for(j in 2:totalProcessors){

          if(( processorsSums[j] + currentValue) < minimum ){

            minimun <- processorsSums[j] + currentValue
            reference <- j

          }

        }}

      # Add the minimum to the list and update everything
      processorsLines[[reference]] <- c(processorsLines[[reference]],analyzeLine)
      processorsSums[reference] <- processorsSums[reference] + currentValue
    }
  }
  # Finally, we have all the lines divided. Create each dataframe and return them
  processorsSubDataframes <- list()
  for(i in 1:totalProcessors){
    #     print(paste("For processor",i))
    #     print(processorsLines[[i]])
    #     print(length(processorsLines[[i]]))
    #     print(processorsSums[i])
    processorsSubDataframes[[i]] <- configDataframe[processorsLines[[i]],]
  }
  return(processorsSubDataframes)
}
