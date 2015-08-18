
.FirstInPair<-function(x){
	intToBits(x)[7]=="01"
}

.SecondInPair<-function(x){
	intToBits(x)[8]=="01"
}

.wrap <- function(...) {
    file.sep <- .Platform$file.sep
    splitted <- paste(unlist(strsplit(paste(...), split = file.sep)),
                                      collapse = paste0(file.sep, " "))
    splitted.pasted <- paste(strwrap(paste(splitted),
                                     exdent = 2), collapse = "\n")
    gsub(paste0(file.sep, " "), file.sep, splitted.pasted)
}
