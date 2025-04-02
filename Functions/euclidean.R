###This function specified a euclidean distance

## Inputs:
## a=vector 1
## b=vector 2

## Outputs:
## euclidean=the euclidean distance between the vectors

euclidean <- function(a, b) sqrt(sum((a - b)^2))