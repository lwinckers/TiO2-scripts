### Function to automatically install packages in R, 
### if they are not yet installed yet (CRAN only).

using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require, character.only=T))
  need <- libs[req==F]
  n <- length(need)
  if(n>0){
    libsmsg <- if(n>2) paste(paste(need[1:(n-1)], collapse=", "), ",", sep="") else need[1]
    if(n>1){
      libsmsg <- paste(libsmsg," and ", need[n], sep="")
    }
    libsmsg <- paste("The following packages could not be found:\r\r",libsmsg,"\r\rInstall missing packages?", sep="")
    if(winDialog(type = c("yesno"), libsmsg)=="YES"){       
      install.packages(need)
      lapply(need, require, character.only=T)
    }
  }
}


