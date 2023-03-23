library(testthat)

# Load data (child to parent relationships)
GRN <- data.frame(rbind(c("bafB,1,bifA"),  c("bifA,1,exbD"), c("expD,1,fur"), c("bupC,2,bafB,bupB"),c("bupB,1,fofH"),c("fofH,1,fur")))
GRN_out <- c("exbD->bifA->bafB->bupC",
             "fur->expD",
             "fur->fofH->bupB->bupC")

#############################################
### Input the text file to generate GRN #####
#############################################

path <- function (GRN){
  relationships <- list()# open the file
  paths <- list()
  # split the line into child and parents
  for (line in 1:nrow(GRN)) {
    line_elements <- strsplit(GRN[line,], ",") [[1]]
     #child <- as.numeric(line_elements[1])
    child <- line_elements[1]
    #child <- GeneIDs[child+1]
    Num.regs <- as.numeric(line_elements[2])
    parents <- strsplit(line_elements[3:(2+Num.regs)]," ")
    #parents <-lapply(parents,as.numeric)
    parents <- unlist(parents)
    #parents <- GeneIDs[parents +1]
    for (parent in parents) {
      if (!(parent %in% names(relationships))) {
        relationships[[parent]] <- list()
      }
      relationships[[parent]] <- c(relationships[[parent]], child)
    }
  }
  #relationships # each child with their parents
  
  # create the paths
  paths <- list()  
  for (parent in names(relationships)) {
    for (child in relationships[[parent]]) {
      path <- child
      current <- child
      while (current %in% names(relationships)) {
        current <- relationships[[current]][[1]]
        path <- paste0(path, "->" ,current)
      }
      path <- paste0(parent,"->",path)
      paths <- c(paths, path)
    }
  }
  
  # Eliminate duplicate paths
  eliminate <- c()
  for (i in 1:length(paths)){
    if(sum(grepl(paths[i], paths[-i], fixed = TRUE))>0 ){
      eliminate <- c(eliminate, i)
    }
  }
  final_paths <- paths[-eliminate]
  final_paths # from children.... to parents...
  
  return (unlist(final_paths))
} 


test_that("single number", {
  expect_equal(path(GRN), GRN_out)
})

