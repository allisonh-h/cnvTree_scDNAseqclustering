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


# loop_tracker() track how many times a target function is called
#===============================================================================
#'
#'
#' @keywords internal
#' 
loop_tracker <- function(f) {
  countLoop <- 0
  function(...) {
    countLoop <<- countLoop + 1
    cat("Function called", countLoop, "times.\n")
    f(...)
  }
}

