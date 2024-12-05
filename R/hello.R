#######################################################################
# Steps for generating Roxygen comments in Rd file
# File -> New project -> Package name (Cancer) -> Create project
# Build (configure Build Tools) -> generate documentation with
#.                                 roxygen ok.
#
# Modify Description (Title, email, name) optional
# To insert roxygen comments skeleton
#.    Code (insert roxygen skeleton)
# Install
# more -> Documents
######################################################

#' Prints "Hello, world!"
#'
#' This function prints the message "Hello, world!" to the console.
#'
#' @return This function has no return value.
#' @export
hello <- function() {
  print("Hello, world!")
}
