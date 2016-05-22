#' This function takes all the files in a given folder that match a regex
#' pattern and put them toguether into a single file.
#'
#' BE CAREFULL: This function __WILL_DELETE__ all the original little files
#' afterward if one of the parameters it set to TRUE.
#'
#' The function takes the following parameters:
#'
#'@param targetFolder (String) The folder where the files are suppose to be.
#'
#'@param regex (String) The files to be combined, must comply with this
#'                           regex.
#'
#'@param finalFileName (String) The name of the final file.
#'
#'@param header (Bool)   If header is true, it will take the header of the
#'                           first file only and ignore the rest. Otherwise it
#'                           will repeat the header many times.
#'
#'@param delete (Bool)   If true, it will delete the original files
#'
#'@param isrecursive (Bool)   If recursive is true, it will look inside all the
#'                           folder of the given folder for more files which
#'                           name match the regex.
#'
#'@return (String) A string with the nucleotides reversed and complemented
unifyFiles <- function (targetFolder, regex, finalFileName, header, delete, isrecursive){

  # Get all the files
  allFiles <- list.files(targetFolder, recursive = isrecursive)

  # Now get those only with the regex expresion
  candidates <- grep(regex, allFiles, fixed=FALSE, ignore.case=FALSE)

  if(length(candidates)>0){
    # Inside this folder it will be a file named as the user desired
    finalFilePath <- paste(targetFolder, "/", finalFileName, sep = '')
    finalFileFD   <- file(finalFilePath, open="w")

    for(i in 1:length(candidates)){
      candidateName <- paste(targetFolder,"/",allFiles[candidates[i]],sep="")
      # Read the file into a string variable
      textFile <- readLines(candidateName,encoding="UTF-8")
      if(header == TRUE && i!=1){
        textFile <- textFile[2:length(textFile)]
      }

      # Write the text into the unified folder
      writeLines(textFile, finalFileFD)

      if(delete == TRUE){
        # Remove the old file
        file.remove(candidateName)
      }
    }
    # Close the file descriptor for the final file
    close(finalFileFD)
  }
  return (length(candidates))
}

#' Delete files from configTables Forward_Reads_File and Reverse_Reads_File (fastq files).
#'
#'@param configTable (data.frame) Contains Forward_Reads_File and Reverse_Reads_File to be removed
#'@return (void) In case of fail, prints err.
deleteFiles <- function(configTable){
  for (i in 1:dim(configTable)[1]){
    file.remove(configTable$Forward_Reads_File[i])
    file.remove(configTable$Reverse_Reads_File[i])
  }
}

#' Unpack fastq files if needed, correct paths to the files in configTable.
#'
#' @param configTable (data.frame) Contains configuration file.
#' @param temp_folder (String) Where to store unzipped files.
#' @return (data.frame) configTable with Forward_File and Reverse_File updated if
#' unpacking files was required.
#' @importFrom R.utils isGzipped gunzip
unpackFastq <- function(configTable, temp_folder){
  for(i in 1:dim(configTable)[1]){
    forward <- configTable$Forward_Reads_File[i]
    rewerse <- configTable$Reverse_Reads_File[i]
    if (isGzipped(forward)) {
      for_name <- sub(".gz", "", forward, fixed=TRUE)
      for_destname <-ifelse(temp_folder != "",
                            paste0(temp_folder, "/", basename(for_name)),
                            for_name)
      gunzip(forward,
             destname = for_destname,
             skip = TRUE, overwrite = FALSE, remove = FALSE)
      configTable$Forward_Reads_File[i] <- for_destname
    }
    if (isGzipped(rewerse)) {
      rew_name <- sub(".gz", "", rewerse, fixed=TRUE)
      rew_destname <- ifelse(temp_folder != "",
                             paste0(temp_folder, "/", basename(rew_name)),
                             rew_name)
      gunzip(rewerse,
             destname = rew_destname,
             skip = TRUE, overwrite = FALSE, remove = FALSE)
      configTable$Reverse_Reads_File[i] <- rew_destname
    }
  }
  return(configTable)
}

#' This function check if the given file exist and can be read.
#'
#' @param filePath (String) A string the path to the file.
#' @return (Void) Stop if no access.
checkFileAccess <- function (filePath){
  if (file.access(filePath[1], mode = 4) != 0) {
    stop(paste0("No read access to the path or it does not exists: ", filePath))
  }
}

#' This function checks if the given directory exist and can be written to.
#'
#' @param filePath (String) A string the path to the file.
#' @return (Void) Stop if no access.
checkFileWriteAccess <- function(filePath){
   if (file.access(filePath[1], mode = 2) != 0) {
     stop(paste0("No write access to the path or it does not exists: ", filePath))
   }
}


#' Gives a config file path and check if exist, is readable, etc...
#'
#'@param configFilePath (string) A string with the absolute path pointing to the config file
#'@return bool, TRUE:  Everything went good.
#'          FALSE: If the config file doesn't exist or if we don't have read access
#'                 If the string provided was 'NULL', 'NA', '0', 'null' or 'na'
checkReady <- function(configFilePath){
  result <- FALSE
  # Check that the config file exist and we have read access
  if(!is.null(configFilePath[1])){ #First we need to check if it is null or we will get a warning
    if(!is.na(configFilePath[1])){ #Same with NA
      # Check that we don't have weird NULL strings
      if(nchar(configFilePath[1]) != 0 && configFilePath[1]!="NULL" && configFilePath[1]!="null" &&
         configFilePath[1] != "na" && configFilePath[1] != "NA"){
        # Check that we have read access, NOTE!: If we do, the function returns 0, and -1 if we don't.
        if(file.access(configFilePath[1], mode = 4) == 0 ){
          result <- TRUE
        }
      }
    }
  }
  return (result)
}
