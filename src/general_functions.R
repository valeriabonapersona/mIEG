## Script: "General functions"

## project: mIEG

## Author: Valeria Bonapersona
## contact: v.bonapersona-2 (at) umcutrecht.nl



# Data handling -----------------------------------------------------------
# Read specific sheets from excel file
read_sheets <- function(filename, sheets, tibble = FALSE) {
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
