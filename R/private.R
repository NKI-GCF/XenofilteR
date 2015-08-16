.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}

## Function to check which reads are first in pair
FirstInPair<-function(x){
	intToBits(x)[7]=="01"
}

## Function to check which reads are second in pair
SecondInPair<-function(x){
	intToBits(x)[8]=="01"
}