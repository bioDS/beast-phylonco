source("simulate_seq.r")

# Simulates a small fixed tree with sequences at the leaves
# and generates a beast runnable xml from the sequence data
# ---------------------------------------------------------------

# simulate tree with newick format
#  (((A:t2, B:t2)D:t1, C:t1 + t2)E);
# lambdaD - substitution model parameter
# lambdaL - substitution model parameter
# t1 - branch times
# t2 - branch times
# l - sequence length
simulateTree <- function(lambdaD, lambdaL, t1, t2, l) {  
  Q <- getQ(lambdaD, lambdaL)
  freq <- getPi(lambdaD, lambdaL)

  e <- generateSeq(l, freq)
  d <- mutateSeq(e, getP(Q, t1))  
  c <- mutateSeq(e, getP(Q, t1 + t2)) 
  b <- mutateSeq(d, getP(Q, t2))
  a <- mutateSeq(d, getP(Q, t2))

  a <- translateSeq(a)
  b <- translateSeq(b)
  c <- translateSeq(c)  

  list(a, b, c)
}

# generate 100 trees with fixed lambdas 
simulateFixed <- function(fileName, lambdaD, lambdaL, seed=777) {
  t1 <- 0.05
  t2 <- 0.05
  l <- 200
  N <- 100

  set.seed(seed)
  dir.create(file.path("../output", "sequences"), showWarnings=F)
  sink(file.path("../output", "sequences", fileName))
  cat("tree,lambdaD,lambdaL,node,sequence\n")
  for (i in 1:N) {
    tree <- simulateTree(lambdaD, lambdaL, t1, t2, l)
    cat(i, lambdaD, lambdaL, "a", "", sep=",")
    cat(tree[[1]], "\n", sep="")
    cat(i, lambdaD, lambdaL, "b", "", sep=",")
    cat(tree[[2]], "\n", sep="")
    cat(i, lambdaD, lambdaL, "c", "", sep=",")
    cat(tree[[3]], "\n", sep="")
  }
  sink()
}

# generate 100 trees with lognormal lambdas
simulateLognorm <- function(fileName, muD, muL, seed=777) {
  t1 <- 0.05
  t2 <- 0.05
  l <- 200
  N <- 100

  set.seed(seed)
  dir.create(file.path("../output", "sequences"), showWarnings=F)
  sink(file.path("../output", "sequences", fileName))
  cat("tree,lambdaD,lambdaL,node,sequence\n")
  for (i in 1:N) {
    lambdaD <- rlnorm(1, muD)
    lambdaL <- rlnorm(1, muL)
    tree <- simulateTree(lambdaD, lambdaL, t1, t2, l)
    cat(i, lambdaD, lambdaL, "a", "", sep=",")
    cat(tree[[1]], "\n", sep="")
    cat(i, lambdaD, lambdaL, "b", "", sep=",")
    cat(tree[[2]], "\n", sep="")
    cat(i, lambdaD, lambdaL, "c", "", sep=",")
    cat(tree[[3]], "\n", sep="")
  }
  sink()
}

# create single xml
createXml <- function(seqData, i, newick, xmlTemplate, dMu, lMu, logName) {
  a <- seqData$sequence[seqData$tree == i & seqData$node == "a"]
  b <- seqData$sequence[seqData$tree == i & seqData$node == "b"]
  c <- seqData$sequence[seqData$tree == i & seqData$node == "c"]
  s <- gsub(pattern="REPLACE_SEQA", replace=a, x=xmlTemplate)
  s <- gsub(pattern="REPLACE_SEQB", replace=b, x=s)
  s <- gsub(pattern="REPLACE_SEQC", replace=c, x=s)
  s <- gsub(pattern="REPLACE_NEWICK", replace=newick, x=s)
  s <- gsub(pattern="REPLACE_LAMBDA_D_MU", replace=dMu, x=s)
  s <- gsub(pattern="REPLACE_LAMBDA_L_MU", replace=lMu, x=s)
  gsub(pattern="REPLACE_NAME", replace=logName, x=s)
}

# create xmls
createXmls <- function(seqName, templateName, newick, dMu, lMu) {
  outputFmt <- sub(".csv", "_NUM.xml", seqName)
  templatePath <- file.path("../templates", templateName)
  seqPath <- file.path("../output", "sequences", seqName)
  xmlTemplate <- readLines(templatePath)  
  seqData <- read.csv(seqPath, colClasses=c("sequence"="character"))
  N <- length(unique(seqData$tree))
  
  dir.create(file.path("../output", "xml"), recursive=T, showWarnings=F)  
  
  for (i in 1:N) {
    fileName <- file.path("../output", "xml", sub("NUM", i, outputFmt))
    logName <- sub(".xml", "", basename(fileName))
    s <- createXml(seqData, i, newick, xmlTemplate, dMu, lMu, logName)
    writeLines(s, fileName)
  }  
}

# run simulations
simulateLognorm("tree_lognorm1.csv", -1, 0.5)
simulateLognorm("tree_lognorm2.csv", 0.5, -1)

# generate xmls
templateName <- "testSiFit3_template.xml"
newick = "((A:0.05, B:0.05)D:0.05, C:0.1)E;"
createXmls("tree_lognorm1.csv", templateName, newick, -1, 0.5)
createXmls("tree_lognorm2.csv", templateName, newick, 0.5, -1)
