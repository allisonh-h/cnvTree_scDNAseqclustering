#' Print Debugging Messages
#'
#' A helper function to print standardized status messages to the console
#' during package execution when the config_hid.yaml `msg` variable is set to TRUE.
#' 
#' @param fucStep A character string indicating the current step or name of the 
#'  function being executed.
#' @param status A character string representing the status (e.g., "start", "end").
#' @param msg A logical value; if \code{TRUE}, the debug message will be printed.
#' 
DebugMsg <- function(fucStep, status, msg)
{
  if (isTRUE(msg)) {
    msg_str <- paste0("LH: ", status, " function", fucStep)
    print(msg_str)
  }
}


#' Record current time for function started timing
#'
#' This utility function records the current system time, typically used
#' for timing the execution of other functions. It is useful for benchmarking
#' or logging the duration of function calls.
#'
#' @param ... Additional arguments passed to methods. Currently not used but included 
#'  for compatibility and extensibility.
#'
#' @return A POSIXct object representing the current system time at the moment this 
#'  function is called.
#'
startTimed <- function(...)
{
  config_path_hid <- system.file("cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" startTimed_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  x <- paste0(..., collapse = "")
  message(x, appendLF = FALSE)
  ptm <- proc.time()
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(ptm)
}


#' Calculate time elapsed since start time
#'
#' This utility function calculates the time elapsed since a recorded start time,
#' typically used for measuring function execution duration. It provides a message
#' displaying the time consumed between two lines of code.
#'
#' @param ptm A POSIXct object representing the start time, typically obtained from 
#'  a call to \code{record_time()}.
#'
#' @return A message showing the time consumed between the start time and the moment 
#'  this function is called.
#'
endTimed <- function(ptm)
{
  config_path_hid <- system.file("cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" endTimed_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
}


#' Export data as a TXT file
#'
#' This function writes a data table to a `.txt` file, saving it to a specified path.
#'
#' @param data A data frame or matrix to be exported as a `.txt` file.
#' @param filename A character string specifying the name of the output `.txt` file.
#' @param path A character string specifying the directory where the file will be saved.
#'
#' @return No return value, called for side effects (writing a file to disk).
#'
writeOutput <- function(data, filename, path)
{
  config_path_hid <- system.file("cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" writeOutput_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  FILEpath <- paste0(path, filename, ".txt")
  utils::write.table(data, file = FILEpath, row.names = FALSE, col.names = TRUE)
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
}
