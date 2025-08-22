utils::globalVariables(c("."))

#' @noRd
seq.split <- function(str, space=50){
  split <- str_split(str, "")[[1]]
  char <- 0
  out.str <- ""
  while(char < length(split)){
    if((char+space)>length(split)){
      end <- length(split)
    }else{
      end <- (char+space)
    }
    if(out.str!=""){
      out.str <- paste0(out.str, "\n")
    }
    out.str <- paste0(out.str, paste(split[(char+1):end], collapse=""))
    char <- char+(space-1)
  }
  return(out.str)
}
