
# From Jeremy
text_to_paths <- function(file_path){
  # create an empty list to store paths
  paths <- list()# create a dictionary to store the parent-child relationships
  relationships <- list()# open the file
  file_content <- readLines(file_path)
  # iterate over the lines in the file
  for (line in file_content) {
    # split the line into child and parents
    line_elements <- strsplit(line, ":")[[1]]
    child <- line_elements[1]
    parents <- strsplit(line_elements[2], " ")[[1]]
    # store the parent-child relationships in the dictionary
    for (parent in parents) {
      if (!(parent %in% names(relationships))) {
        relationships[[parent]] <- list()
      }
      relationships[[parent]] <- c(relationships[[parent]], child)
    }
  }
  # create the paths
  for (parent in names(relationships)) {
    for (child in relationships[[parent]]) {
      path <- child
      current <- child
      while (current %in% names(relationships)) {
        current <- relationships[[current]][[1]]
        path <- paste0(current,path)
        paths <- c(paths, path)
      }
    }
  }
}

file_path = "graph.txt"
paths <- text_to_paths(file_path)

###############
### My code ###
###############

# Load data
GRN <- read.table("~/Missing data/SERGIO/Demo/single_cell_GRN.txt", quote="\"", comment.char="")
GRN_out <- read.csv("~/Missing data/SERGIO/Demo/GeneExpression.csv")

#GRN <- data.frame(rbind(c("bafB,1,bifA"),  c("bifA,1,exbD"), c("expD,1,fur"), c("bupC,2,bafB,bupB"),c("bupB,1,fofH"),c("fofH,1,fur")))

###############
### Input #####
###############

relationships <- list()# open the file
paths <- list()
# split the line into child and parents
for (line in 1:nrow(GRN)) {
  line_elements <- strsplit(GRN[line,], ",") [[1]]
  #child <- as.character(as.numeric(line_elements[1]))
  child <- line_elements[1]
  child <- paste0("Gene",child)
  Num.regs <- as.numeric(line_elements[2])
  parents <- strsplit(line_elements[3:(2+Num.regs)]," ")
  #parents <-lapply(parents,as.numeric)
  parents <-lapply(parents,as.character)
  for (parent in parents) {
    if (!(parent %in% names(relationships))) {
      relationships[[parent]] <- list()
    }
    relationships[[parent]] <- c(relationships[[parent]], child)
    #relationships[[parent]] <-lapply(relationships[[parent]],as.numeric)
  }
}
relationships

# create the paths
paths <- list()  
for (parent in names(relationships)) {
  for (child in relationships[[parent]]) {
    path <- child
    current <- child
    while (current %in% names(relationships)) {
      print(current[[1]])
      current <- relationships[[current]][[1]]
      print(current[[1]])
      path <- paste0(current,"-" ,path)
      print("itr")
    }
    path <- paste0(path,"-",parent)
    print("itrchild")
    paths <- c(paths, path)
  }
}
paths

# Eliminate duplicate paths
eliminate <- c()
for (i in 1:length(paths)){
  if(sum(grepl(paths[i], paths[-i], fixed = TRUE))>0 ){
    eliminate <- c(eliminate, i)
  }
}
final_paths <- paths[-eliminate]
final_paths # from children.... to parents...


###############
### Output ####
###############

summary(GRN_out)
dim(GRN_out)
res <- cor(t(GRN_out[,-1]))
corrplot::corrplot(res,type = "upper", tl.col = "black", tl.cex = 0.5)


