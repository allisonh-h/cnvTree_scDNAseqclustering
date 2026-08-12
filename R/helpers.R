#' Initialize Package Options on Load
#'
#' Loads internal configuration from the cnvTree_config_hid.yanl and sets global 
#' options.
#'
#' @param libname Library name (passed automatically by R)
#' @param pkgname Package name (passed automatically by R)
#'
#' @noRd
#' 
.onLoad <- function(libname, pkgname) 
{
  config_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = pkgname)
  cnvTree_config_hid <- yaml::read_yaml(config_path)
  
  # store the raw values from the cnvTree_config_hid.yaml file
  options(cnvTree_v_num      = cnvTree_config_hid$cnvTree_v_num)
  options(cnvTree_msg        = cnvTree_config_hid$cnvTree_msg)
  options(sexchromosome      = cnvTree_config_hid$sexchromosome)
  options(cluster            = cnvTree_config_hid$cluster)
  options(smoothheatmap      = cnvTree_config_hid$smoothheatmap)
  options(min_cell           = cnvTree_config_hid$min_cell)
  options(overlap_region     = cnvTree_config_hid$overlap_region)
  options(dif_ratio          = cnvTree_config_hid$dif_ratio)
  options(consecutive_region = cnvTree_config_hid$consecutive_region)
  options(difratio_chr       = cnvTree_config_hid$difratio_chr)
  options(min_cell_subclone  = cnvTree_config_hid$min_cell_subclone)
}


#' DebugMsg
#'
#' Print Debugging Messages; a helper function to print standardized status messages 
#' to the console during package execution when the config_hid.yaml `msg` variable 
#' is set to TRUE.
#' 
#' @param fucStep A character string indicating the current step or name of the 
#'  function being executed.
#' @param status A character string representing the status (e.g., "start", "end").
#' @param cnvTree_msg A logical value; if \code{TRUE}, the debug message will be printed.
#' 
#' @noRd
#' 
DebugMsg <- function(fucStep, status, cnvTree_msg)
{
  if (isTRUE(cnvTree_msg)) {
    msg_str <- paste0("LH: ", status, " function", fucStep)
    print(msg_str)
  }
}


#' startTimed
#'
#' Record current time for function started timing; this utility function records 
#' the current system time, typically used for timing the execution of other functions. 
#' It is useful for benchmarking or logging the duration of function calls.
#'
#' @param ... Additional arguments passed to methods. Currently not used but included 
#'  for compatibility and extensibility.
#'
#' @return A POSIXct object representing the current system time at the moment this 
#'  function is called.
#'
#' @noRd
#' 
startTimed <- function(...)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" startTimed_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  x <- paste0(..., collapse = "")
  message(x, appendLF = FALSE)
  ptm <- proc.time()
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(ptm)
}


#' endTimed
#'
#' Calculate time elapsed since start time; this utility function calculates the 
#' time elapsed since a recorded start time, typically used for measuring function 
#' execution duration. It provides a message displaying the time consumed between 
#' two lines of code.
#'
#' @param ptm A POSIXct object representing the start time, typically obtained from 
#'  a call to \code{record_time()}.
#'
#' @return A message showing the time consumed between the start time and the moment 
#'  this function is called.
#'  
#' @noRd
#' 
endTimed <- function(ptm)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" endTimed_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}


#' writeOutput
#'
#' Export data as a TXT file; this function writes a data table to a `.txt` file, 
#' saving it to a specified path.
#'
#' @param data A data frame or matrix to be exported as a `.txt` file.
#' @param filename A character string specifying the name of the output `.txt` file.
#' @param path A character string specifying the directory where the file will be 
#'  saved.
#'
#' @return No return value, called for side effects (writing a file to disk).
#'
writeOutput <- function(data, filename, path)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" writeOutput_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------

  FILEpath <- paste0(path, filename, ".txt")
  
  #if (is.list(data)) { #LH added 072026 -->>
  #  data <- sapply(data, function(x) paste(unlist(x), collapse = ", "))
  #}
  #LH added 072026 <<--
  utils::write.table(data, file = FILEpath, row.names = FALSE, col.names = TRUE)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}
