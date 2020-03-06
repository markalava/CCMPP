#' ccmpp: Cohort component method of population projection.
#'
#' @docType package
#' @name ccmpp
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("ccmpp", libpath)
}
