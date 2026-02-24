# DebugMsg() Print debug messages for internal tracking
#'==============================================================================
#' A helper function to output stage-specific messages to the console 
#' when the global `msg` variable is set to TRUE.
#' 
#' @param fucStep A character string or numeric indicating the specific 
#' step within the function.
#' @param stage A character string indicating the current processing 
#' phase (e.g., "Preprocessing", "Clustering").
#' 
#' @return None. Prints a formatted string to the console.
#' 
#' @details 
#' The output is prefixed with "LH: " to easily identify messages 
#' originating from this package during a noisy console session.
#' 
#' @keywords internal
#'
DebugMsg <- function(fucStep, stage)
{
  if (msg) {
    msg_str <- paste0("LH: ", stage, " function", fucStep)
    print(msg_str)
  }
}


# startTimed(), endTimed() for time calculating
#===============================================================================
#' Record current time for function started timing
#'
#' This utility function records the current system time, typically used
#' for timing the execution of other functions. It is useful for benchmarking
#' or logging the duration of function calls.
#'
#' @param ... Additional arguments passed to methods.
#' Currently not used but included for compatibility and extensibility.
#'
#' @return A POSIXct object representing the current system time at the moment this function is called.
#'
#' @keywords internal
#'
startTimed <- function(...){
  fucStep <- " utils_startTimed_cnvTree_v030"
  DebugMsg(fucStep, "start")
  x <- paste0(..., collapse = "")
  message(x, appendLF = FALSE)
  ptm <- proc.time()
  DebugMsg(fucStep, "end")
  return(ptm)
}


# endTimed: Calculate time elapsed since start time
#'
#' This utility function calculates the time elapsed since a recorded start time,
#' typically used for measuring function execution duration. It provides a message
#' displaying the time consumed between two lines of code.
#'
#' @param ptm A POSIXct object representing the start time, typically obtained from a call to \code{record_time()}.
#'
#' @return A message showing the time consumed between the start time and the moment this function is called.
#'
#' @keywords internal
#'
endTimed <- function(ptm){
  fucStep <- " utils_endTimed_cnvTree_v030"
  DebugMsg(fucStep, "start")
  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")
  DebugMsg(fucStep, "end")
}


# writeOutput() function for outputting cluster results in .txt
#'==============================================================================
#' Export data as a TXT file
#'
#' This function writes a data table to a `.txt` file, saving it to a specified path.
#'
#' @param data A data frame or matrix to be exported as a `.txt` file.
#' @param filename A character string specifying the name of the output `.txt` file.
#' @param path A character string specifying the directory where the file will be saved.
#'
#' @return The function writes a `.txt` file and returns the file path as a character string.
#' @export
#'
writeOutput <- function(data, filename, path){
  fucStep <- " utils_writeOutput_cnvTree_v030"
  DebugMsg(fucStep, "start")
  FILEpath <- paste0(path, filename, ".txt")
  utils::write.table(data, file = FILEpath, row.names = FALSE, col.names = TRUE)
  DebugMsg(fucStep, "end")
}
