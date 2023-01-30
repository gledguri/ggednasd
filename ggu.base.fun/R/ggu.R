#' Generates colours that are not white
#'
#' This function randomly generates 'n' colours that are not white
#' @param n The number of colours to be produced
#' @return R vector
#' @export
#' @examples
#' colur(4)
#' [1] "darkorange2" "orchid1"     "skyblue1"    "mediumblue"

colur<-function(n){
  temp<-as.vector(colors())
  temp1<-temp[sample(c(8:151,362:657),n, replace=FALSE)]
  return(temp1)
}

#' Backtransforms the logit values
#'
#' This function transforms values from logit into values between 0 and 1
#' @param n The number/vector that has logit values
#' @return R number or vector with backtransformed logit values
#' @export
#' @examples
#' inverselogit(4)

inverselogit <- function(x){
  return(1/(1+exp(-x)))
}

#' Transforms values into logit-scale
#'
#' This function transforms values between 0 and 1 into logit-scale values
#' @param n The number/vector that has [0:1] values
#' @return R number or vector with logit values
#' @export
#' @examples
#' logit(4)
logit <- function(x){
  return(log(x/(1-x)))
}


#' Generates string text with colours
#'
#' This function returns a text with colourw
#' @param input the text string to be printed
#' @param col the code for the color
#' @return text string
#' @export
#' @examples
#' questions("hehehe")
#' hehehe (colored)


questions <- function(input,col=43){
  cat(paste0("\033[0;", col, "m", input,"\033[0m"))
}


#' Prints out the statistical test in a .txt file
#'
#' This function prints the statistical results in a text file with the name of the object as the name of the file.
#' @param title The name of the R object
#' @param wd wd=getwd(), if you want to change this you could write paste0(wd,"/",folder)
#' @return Creates a .txt file in the wd(working directory) with the results
#' @export
#' @examples
#' write.result(t.test(mydata))
#' write.result(adonis(dist ~ groups, permutations= 999))

write.result<-function(title, wd=getwd()){
  x<-paste0(wd,"/",substitute(title),".txt")
  cat("\n--------------------------------------------------------------------------\n", file = x, append = TRUE)
  cat(substitute(title), file = x)
  cat("\n--------------------------------------------------------------------------\n", file = x, append = TRUE)
  strwrap(capture.output(title, file = x, append = TRUE))
  return(x)
}


#' Converts a species abundance table into relative abundance (*100%)
#'
#' This function converts the abundances into relative abundance per station (column)
#' @param input Species abundances R object (species are rows x stations are columns)
#' @param threshold threshold=0, if you want to change the threshold where the % will be considered 0
#' @return Creates a R object with the same column & row names as the input
#' @export
#' @examples
#' output<-rel.ab(sp.abundance.data) #columns are stations and rows are species
#' output<-rel.ab(sp.abundance.data, threshold = 0.5)

rel.ab<-function(input,threshold=0){
  a<-data.frame(matrix(NA,nrow(input),ncol(input)))
  colnames(a)<-colnames(input)

  for (i in 1:nrow(input)){
    a[i,]<-input[i,]/colSums(input)*100
  }
  rownames(a)<-rownames(input)
  a[a<threshold]<-0
  return(a)
}



#' Collapse otu table based on chosen level
#'
#' This function for collapsing otu statoins (columns) based on n=column name in metadata
#' @param otu otu table where species are rows and columns are stations
#' @param metadata metadata where all stations are rows
#' @return collapsed otu table based on chosen level
#' @export
#' @examples
#' otu<-collapse_otu #columns are stations and rows are species

collapse_otu<-function(otu,metadata=mmn){
  cat("Is collapse.name column created in the metadata?")
  bb<-as.data.frame(otu)
  colnames(bb) <- metadata[,"collapse.name"]
  otu_coll <- data.frame(matrix(0,nrow(bb),length(unique(colnames(bb)))))
  colnames(otu_coll) <- unique(colnames(bb))
  for (i in unique(colnames(bb))) {
    otu_coll[,colnames(otu_coll)==i] <- rowSums(bb[colnames(bb)==i])
  }
  return(as.data.frame(otu_coll))
}


#' Collapse metadata table based on chosen level
#'
#' This function for collapsing otu statoins (columns) based on n=column name in metadata
#' @param metadata metadata where all stations are rows
#' @return collapsed metadata table based on chosen level
#' @export
#' @examples
#' mm<-collapse_meta #remember to collapse on the column "collapse.name"

collapse_meta<-function(metadata){
  metadata[!duplicated(metadata[,"collapse.name"]),]
}


#' Collapse metadata table based on chosen level
#'
#' This function for collapsing otu statoins (columns) based on n=column name in metadata
#' @param otu otu table where species are rows and columns are stations
#' @param threshold the threshold where a MOTU is considered "Present"
#' @return otu table with only presence / absence
#' @export
#' @examples
#' otu<-pa_conversion(otu)
pa_conversion <- function(otu,threshold=0){
  threshold <- threshold
  for (i in 1:ncol(otu)) {
    for (j in 1:nrow(otu)) {
      if(otu[j,i]>threshold){
        otu[j,i] <- 1}
      else if(otu[j,i]<=threshold){otu[j,i]<-0}
    }
  }
  return(otu)
}
#' Creates a piechart on MOTU diversity
#'
#' This function creates a pie-chart on the diversity of MOTUs found grouped by the "rank"
#' @param taxa Dataframe that contains taxa (each rank in one column)
#' @param rank The "rank" that you are interested in grouping the MOTUs to
#' @param logic If logic=2 the "unassigned" MOTUs will be removed. If logic=1 the "unassigned" MOTUs will remain
#' @return Creates a pie-chart with the diversity of MOTUs per "rank"
#' @export
#' @examples
#' motudiv(taxa,"species")
#' motudiv(taxa,"phylum_name")
motudiv<-function(taxa, rank, logic=1){

require(plyr,dependencies = TRUE)
require(stringr,dependencies = TRUE)
require(ggplot2,dependencies = TRUE)
  
  xx<-logic
  x<-as.data.frame(taxa[,rank])
  otu.pie<-data.frame(x,c(1))
  colnames(otu.pie)<-c(rank, "reads")
  otu.pie<-ddply(otu.pie, c(rank), summarize, sp=sum(reads))
  if(xx==2){
    otu.pie[,1]<-as.character(otu.pie[,1])
    a<-str_detect(otu.pie[,1], "")
    otu.pie<-otu.pie[a,]
  }
  else if(xx==1){
    otu.pie[,1]<-as.character(otu.pie[,1])
    otu.pie[!str_detect(otu.pie[,1], ""),1]<-c("Unassigned")
    otu.pie[,1]<-as.factor(otu.pie[,1])
  }
  otu.pie[,1]<-paste0(otu.pie[,1]," ",otu.pie[,2])
  otu.pie$perc<-round(otu.pie[,2]/sum(otu.pie[,2]),3)
  print(otu.pie)

  ggplot(otu.pie, aes(x = "", y=sp, fill = otu.pie[,rank])) +
    #geom_bar(width = 1, stat="identity")+
    geom_col(width=1)+
    coord_polar(theta = "y")+
    geom_bar(stat="identity", width=2, color="white") +
    theme_void()+
    ggtitle(paste0("The number of different MOTUs detected grouped by ",rank))
}


#' The negative of %in%
#'
#' The negative of %in%
#' @return a vector
#' @export
#' @examples
#' "he"%notin%"hehehe"

`%notin%` <- Negate(`%in%`)


#' Set up a working directory (easily to convert between windows and mac)
#'
#' This function creates a pie-chart on the diversity of MOTUs found grouped by the "rank"
#' @param comp Either "mac" or "win" for specifying in which computer are you working atm
#' @param folder The path to the folder (After the path to OneDrive) that you want to set the working directory
#' @return generates a working directory
#' @export
#' @examples
#' work("IMR/Data/Trawl data")
#' work("IMR/Data")

work<-function(folder=""){
  if(Sys.info()["nodename"]=="HI-13521.hi.no"){
    pd <- (paste0("/Users/","a36142","/OneDrive - Havforskningsinstituttet",folder))
    return(pd)
  }
  else if(Sys.info()["nodename"]=="GLED-WIN"){
    pd <- (paste0("/Users/","gledg","/OneDrive - Havforskningsinstituttet",folder))
  return(pd)
  }
}



#' Set up a working directory (easily to convert between windows and mac)
#'
#' This function creates a pie-chart on the diversity of MOTUs found grouped by the "rank"
#' @param folder The path to the folder (After the path to pd) that you want to set the working directory
#' @param m logical value if you want the menu to list all the directories below that folder
#' @return generates a working directory
#' @export
#' @examples
#' low(folder="eDNA",m=F)
#' low(folder="",m=T)
#'
low <- function(folder="",m=F){
  if(m==T){
    cat(questions("Where do you want to set the working directory"))
    sl <- menu(list.dirs(folder, recursive = F))
    return(paste0((list.dirs(folder, recursive = F)[sl])))
  }else{
    return(paste0(pd,"/",folder))
  }
}



#' Remove blank rows from OTU and TAX
#'
#' This function removes the blanks from OTU and TAX dataframe
#' @param df the input dataframe which the rows will be removed
#' @param c the logical vector expressing which rows should be removed
#' @return Removes the the c rows
#' @export
#' @examples
#' tax<-rr(tax,!c) #removing the c rows from tax data.frame
#' otu<-rr(otu,!c) #removing the c rows from otu data.frame

rr<-function(df,c){
  df<-df[c,]
  return(df)
}


#' Change the parameters in the plot view
#'
#' Change the parameters in the plot view
#' @param x,y input x and y
#' @return the plot window is separated in x rows and y columns
#' @export
#' @examples
#' v(2,2)

v <- function(x=1,y=1){
  par(mfrow=c(x,y))
}


#' Transposes a data.frame from 2 dimensional to 1 dimensional
#'
#' This function creates a pie-chart on the diversity of MOTUs found grouped by the "rank"
#' @param df The data.frame that has the vlaues in the cell as a combination of col and rows
#' @return generates a data.frane with $X = rows, $Y = columns and $Z = the value of the cells
#' @export
#' @examples
#' trans(edna)
#' trans(abundance.data)

trans<-function(df){
  x<-rownames(df)
  y<-colnames(df)
  data <- expand.grid(X=x, Y=y)
  data$Y<-as.character(data$Y)
  data$Z<-NA
  for (i in 1:length(y)){
    n<-((i-1)*length(x))+1
    j<-n+length(x)-1
    co<-y[i]
    data$Z[n:j]<-unlist(df[co])
  }
  return(data)
}


#' Copies the rownames and columnames of a data.frame into another one
#'
#' This function creates a pie-chart on the diversity of MOTUs found grouped by the "rank"
#' @param df The data.frame that the values will be retreived
#' @return generates a data.frane with the same colnames and rownames
#' @export
#' @examples
#' copy(edna)
#' copy(taxa)

copy<-function(df){
  df2<-as.data.frame(matrix(NA,nrow = nrow(df), ncol = ncol(df)))
  colnames(df2)<-colnames(df)
  rownames(df2)<-rownames(df)
  return(df2)
}



#' Joins values from another dataframe
#'
#' This function matches 2 vectors of different dataframes and joins the values attached to the second vector
#' @param input.x The vector which the values will be joined to
#' @param input.y The bridge vector. The vector that has the same as the first vector and rows correspond to values that will be joined
#' @param input.yy The vector which contains the values that will be joined
#' @return generates a vector with only joined values ordered for corresponding to the first vectgor
#' @export
#' @examples
#' mm$category <- merr(mm$samp_cat,category.transformation$samp_cat,category.transformation$category)

merr <- function(input.x,input.y,input.yy){
  yy <- vector(mode="character", length = length(input.x))
  for (i in 1:length(input.x)) {
    yy[i] <- input.yy[input.x[i]==input.y]
  }
  return(yy)
}





#' Print bold
#'
#' This function prints on Console with ##
#' @param input the text that will be printed
#' @param wide How wide should be the print result
#' @param col the colur of the text in between ##
#' @return prints on Console
#' @export
#' @examples
#' printbold("hehehe")
#' printbold("hehehe", wide=60)
#' printbold("hehehe", wide=60, col=38)

printbold <- function(input,wide=40,col=47){cat(cat("\n"),cat(rep("#",wide),sep = ""),
                                                cat("\n"),cat(rep("#",((wide-2) - (nchar(input)))/2),sep=""),cat(" "),
                                                cat(paste0("\033[0;", col, "m", input,"\033[0m")),cat(" "),cat(rep("#",((wide-2) - (nchar(input)))/2),sep=""),
                                                cat("\n"),cat(rep("#",wide),sep = ""))
}


#' Text colors
#'
#' This function prints all text colors
#' @return prints on all text colors
#' @export
#' @examples
#' hh.textcolor()

hh.textcolor <- function()
for (i in 29:49) {
  cat(paste0(questions("colours", i)," - ",i));cat("\n")
}


#' Inverse logit function
#'
#' This function returns the inverse logit value
#' @param x The value in logit scale
#' @return The vale on the normal scale
#' @export
#' @examples
#' inverselogit(0.01)

inverselogit <- function(x){
  return(1/(1+exp(-x)))
}

#' Logit function
#'
#' This function returns the logit value
#' @param x The value in normal scale
#' @return The vale on the logit scale
#' @export
#' @examples
#' logit(0.01)
logit <- function(x){
  return(log(x/(1-x)))
}

#' Match multiple inputs (grep in for multiple variables at once)
#'
#' This function serches and matches multiple inputs at once
#' @param pattern The pattern that that you want to search for. Can be a single variable or a vector
#' @param x The variable or the vector that you want to look in
#' @return The matching values
#' @export
#' @examples
#' match_multiple("Squalius acanthius",dataframe)
#' match_multiple(c("sp1","sp2"),vector)
match_multiple <- function(pattern, x){
  unique (grep(paste(pattern,collapse="|"), x, value=TRUE))
}


#' Multigrep
#'
#' This function does grep but allows for multiple inputs (so a vector)
#' @param pattern The pattern that that you want to search for. Can be a single variable or a vector
#' @param x The variable or the vector that you want to look in
#' @return The matching values
#' @export
#' @examples
#' match_multiple("Squalius acanthius",dataframe)
#' match_multiple(c("sp1","sp2"),vector)
multi_grep <- function(pattern, x){
  grep(paste(pattern,collapse="|"), x)
}

#' Print as dataframe
#'
#' This function prints a multiple vectors in the form of a dataframe. Mainly for esthetic reasons
#' @param vector_var The variable/vector that you want to print
#' @param wide The number of wide space you want from the first variable to the next 
#' @return The formated output in Console
#' @export
#' @examples
#' pr_df(tax$id,50)
#' pr_df(tax$sequences,50,new_row = T)
#' sliced_sequences <- pr_df(tax$sequence,50) #if you want to slice all the variables in that vector to 50 characters 
#' cat(pr_df(tax$genus_name,30),pr_df(tax$species_name,30),pr_df(tax$taxid,10),sep = "\n")
pr_df <- function(vector_var,wide=40,new_row=F){
  temp <- wide-nchar(vector_var)
  
  if (is.null(dim(vector_var))) {
  for (i in 1:length(temp)) {
    if (temp[i]>=0) {
      temp[i] <- paste(rep(" ",temp[i]), collapse = "")
    }else{
      vector_var[i] <- substr(vector_var[i],0,wide)
      temp[i] <- ""
    }
  }
  if (new_row==T) return(cat(paste0(vector_var,temp),sep = "\n"))
  return(paste0(vector_var,temp))
  }
}

#' Clean OTU table from empty rows
#'
#' This function removes the empty rows from the OTU table and is also removes them from 
#' the TAX dataframe. This functions automatically updates the otu and tax table as it creates
#' 2 objects named otu and tax. 
#' It will print out which tax was removed
#' @param otu the input dataframe (eDNA reads) which the rows will be removed
#' @param tax the input dataframe (taxonomic information) which the rows will be removed
#' @return OTU table without rows that contains 0's across all stations
#' @export
#' @examples
#' clean.otu(otu,tax)
#' taxa_removed <- clean.otu(otu,tax)

clean.otu<-function(otu,tax){
  tax<<-tax[rowSums(otu)>0,]
  otu<<-otu[rowSums(otu)>0,]
}

#' Inverse transpose dataframe
#'
#' This function does the inverse of trans() function. It transposes 1 dimentional dataframe
#' into a 2 dimentional dataframe 
#' @param input.X the vector where the rows/rownames are stored
#' @param input.X the vector where the columns/colnames are stored
#' @param input.Z the vector where the values are stored
#' @return A two dimentional dataframe
#' @export
#' @examples
#' clean.otu(otu,tax)
#' taxa_removed <- clean.otu(otu,tax)
inv.trans <- function(input.X,input.Y,input.Z){
  output <- as.data.frame(matrix(NA,length(unique(input.X)),length(unique(input.Y))))
  rownames(output) <- unique(input.X)
  colnames(output) <- unique(input.Y)
  for (i in rownames(output)) {
    for (j in colnames(output)) {
      output[i,j] <- input.Z[input.X==i&input.Y==j]
    }
  }
  return(output)
}

#' Sequence developer
#'
#' This function creates a sequence of numbers from 1:k when i=0 and k+1:2k 
#' This is very handy when you want to slice data into fragments 
#' @param k the length of the sequence of numbers
#' @param i the loop factor
#' @return A vector of a sequence of numbers k-long
#' @export
#' @examples
#' seq_dev(100,0)
#' seq_dev(100,1)
#' for(i in 0:3) vector[i] <- seq_dev(100,i)
seq_dev <- function(k,i) ((k*i)+1):((k*i)+k)

#' eDNA data to fasta
#'
#' This function creates fasta sequences from edna reads
#' This is very handy when you want to reblast some sequences
#' @param qid The vector containig the query ID
#' @param seq The vector containig the sequences
#' @param output the name of the output file
#' @return A fasta file with extenxtion .fasta
#' @export
#' @examples
#' # edna_to_fasta(edna$id,edna$sequence,"Plate_NANB.fasta")
if (!require("seqinr")) {install.packages("seqinr",dependencies = TRUE);require("seqinr")}
edna_to_fasta <- function(qid,seq,output="out.fasta"){
  capture.output(
    for (i in 1:length(qid)) writeLines(paste0(">query",qid[i],";","\n",seq[i])),
    append = F,type = "output",file = "temp.fasta")
  st <- read.fasta("temp.fasta")
  write.fasta(st, qid,file.out = output);unlink("temp.fasta")
}

#' Row Paste
#'
#' This function pastes rows one by one
#' @param inn the dataframe where you want the columns should be collapsed
#' @return A vector of pasted materials
#' @export
#' @examples
#' rowpaste(tax[,1:6])
#' metazoans <- rowpaste(tax[,tax$superkindom=="Metazoa"])
rowpaste <- function(inn){
  vvv <- vector(length=length(nrow(inn)))
  for (k in 1:nrow(inn)) {
    vvv[k] <- paste(inn[k,],collapse = " ")
  }
  return(vvv)
}


#' Annotate blast results
#'
#' This function reads blasts results and automatically or manually annotates the blasted sequences
#' to the best match resutl
#' @param input the blast output dataframe from NCBI
#' @param method The annotation method. 1=Fully manual; 2=Semi-automatic; 3=Fully Automatic
#' @return prints a csv file in your working directory named "manually_annotated_function_outcome.csv"
#' @return returns a dataframe containing query id, species names of query match and % match
#' @export
#' @examples
#' manually_annotate(blast_output,method = 3)
#' @name Packages asd

manually_annotate <- function(input,method=1,skip=T,print_blast=T,query_vector=NA){ #Method 1: "Conservative; Methood 2: "Semiconservative"; Method 3: "Automatic"

if (!require("devtools")) {install.packages("devtools",dependencies = TRUE);require("devtools")}
if (!require("dplyr")) {install.packages("dplyr",dependencies = TRUE);require("dplyr")}
if (!require("coda")) {install.packages("coda",dependencies = TRUE);require("coda")}
if (!require("vctrs")) {install.packages("vctrs",dependencies = TRUE);require("vctrs")}
if (!require("stringr")) {install.packages("stringr",dependencies = TRUE);require("stringr")}

  magic.pident <- function(input,pident){
    input[pident>=99.5] <-             paste0("\033[0;", 32, "m", input[pident>=99.5],"\033[0m")
    input[pident<99.5&pident>=98.5] <- paste0("\033[0;", 36, "m", input[pident<99.5&pident>=98.5],"\033[0m")
    input[pident<98.5&pident>=97.5] <- paste0("\033[0;", 34, "m", input[pident<98.5&pident>=97.5],"\033[0m")
    input[pident<97.5&pident>=96.5] <- paste0("\033[0;", 33, "m", input[pident<97.5&pident>=96.5],"\033[0m")
    input[pident<96.5] <-             paste0("\033[0;", 31, "m", input[pident<96.5],"\033[0m")
    return(input)}
  
  function_col <- function(x = character()) {
    vec_assert(x, character())
    new_vctr(x, class = "vctrs_function_col")}
  
  format.vctrs_function_col <- function(x,...) {
    gsub("function",crayon::red("function"),vec_data(x))}
  
  pr_df <- function(vector_var,wide=40){
    wide=wide
    temp <- wide-nchar(vector_var)
    for (i in 1:length(temp)) {
      if (temp[i]>=0) {
        temp[i] <- paste(rep(" ",temp[i]), collapse = "")
      }else{
        vector_var[i] <- substr(vector_var[i],0,wide)
        temp[i] <- ""
      }
    }
    return(paste0(vector_var,temp))
  }
  
  allduplicates <- function(vector){
    vector%in%unique(vector[duplicated(vector)])
  }
  
  unlist_uneven_list <- function(list){
    v <- vector("numeric")
    for (i in 1:length(list)) {
      v[i] <- length(list[[i]])
    }
    x <- as.data.frame(matrix(NA,length(list),max(v)))
    for (i in 1:length(list)) {
      x[i,1:length(list[[i]])] <- list[[i]]
    }
    return(x)
  }
  
  is.defined <- function(sym) {
    sym <- deparse(substitute(sym))
    env <- parent.frame()
    exists(sym, env)
  }
  
  rowpaste <- function(inn){
    vvv <- vector(length=length(nrow(inn)))
    for (k in 1:nrow(inn)) {
      vvv[k] <- paste(inn[k,],collapse = " ")
    }
    return(vvv)
  }
  
  input <- tibble::tibble(input)
  columns <- c("Ssciname","scommname","qseqid","sseqid","pident",
               "length","mismatch","gapopen","qcovus","qstart","qend",
               "sstart","send","evalue","staxids","qlen","qcovs")
  if (!"annotated_tax"%in%colnames(input)) input$annotated_tax <- NA
  if (!"pmatchsel"%in%colnames(input)) input$pmatchsel <- NA
  cat("The data frame should contain these columns \n",columns)
  cat("checking the names of the input column")
  col_missing <- sum(!columns%in%colnames(input))
  if(col_missing>0)cat(columns[!columns%in%colnames(input)]," column is missing")
  input$Ssciname[grep(";",input$Ssciname)] <- gsub(";", " | ", input$Ssciname[grep(";",input$Ssciname)])
  if(!sum(is.na(query_vector))){input <- input[input$qseqid%in%query_vector,]}
  
  if (method==1) {
    j=1 #build the loop
    while (j <= length(unique(input$qseqid))){ #Loop of method 1
      pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
      i <- unique(input$qseqid)[j] #select the unique query ID
      list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
      
      if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 1
      else { #Not skipping method 1
        #Esthetics of the list
        list_of_query_print <- tibble::tibble(list_of_query);
        list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
        list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
        list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
        list_of_query_print$pident <- function_col(list_of_query_print$pident)
        
        #Print
        if (print_blast==T) {
          cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
          print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
        }
        
        #Create a summary list
        summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
        summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
        
        #Menu and choices
        chsl <- c(paste0(summary_of_list$Ssciname),"Undefined","Enter manually") #Create taxa selection options
        cat("\n");questions("Select the subject that your query belongs to ::");cat("\n")
        sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                                 magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                                 magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
                           "Unidentified","Enter manually"))
        
        if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa ")}
        #Writing on the fasta file
        input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") #Assigns the selected taxa to the input database
        input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
        
        #Esthetics
        setTxtProgressBar(pb,j);cat("\n") #progress bar
        j <- j+1 #loop mechanism
        write.csv(input,"manually_annotated_function_outcome.csv",row.names = FALSE)
      } #Not skipping method 1
    } #Loop of method 1
  } #Method 1
  else if (method==2) { #Method 2
    j=1 #build the loop
    while (j <= length(unique(input$qseqid))){ #Loop of method 2
      pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
      i <- unique(input$qseqid)[j] #select the unique query ID
      list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
      
      #Esthetics of the list
      if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 2
      else { #Not skipping method 2
        list_of_query_print <- tibble::tibble(list_of_query);
        list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
        list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
        list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
        list_of_query_print$pident <- function_col(list_of_query_print$pident)
        
        #Print
        if (print_blast==T) {
          cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
          print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
        }
        
        #Create a summary list
        summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
        summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
        
        chsl <- c(paste0(summary_of_list$Ssciname),"Undefined","Enter manually") #Create taxa selection options
        
        match_id <- 100
        l <- summary_of_list$pident_max==match_id
        if (sum(l)==0) {
          l <- summary_of_list$pident_max==max(summary_of_list$pident_max)
          match_id <- max(summary_of_list$pident_max)
        }
        if (sum(l)==1) {#If only 1 subject is 100%
          cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
          cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                      magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                      magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
          cat("\n");questions("The algorithm selected ::")
          cat("\n");cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname[l],30),summary_of_list$pident_max[l]),"::",
                                magic.pident(pr_df(summary_of_list$pident_max[l],5),summary_of_list$pident_max[l]), 
                                magic.pident(pr_df(summary_of_list$pident[l],5),summary_of_list$pident[l]))),sep="\n")
          cat("\n");questions("Do you agree ?")
          slyn <- menu(c("yes","no"))
          if (slyn==1) {#If yes agree to algorithm's suggestion on annotating to the only subject that is 100#
            input$annotated_tax[input$qseqid==i] <- paste(summary_of_list$Ssciname[l])
            input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%summary_of_list$Ssciname[l]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
            
          }#If yes agree to algorithm's suggestion on annotating to the only subject that is 100#
          else if (slyn==2) {#If no agree to algorithm's suggestion on annotating to the only subject that is 100#
            cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
            sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                                     magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                                     magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
                               "Unidentified","Enter manually"))
            
            if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa ")}
            
            #Writing on the fasta file
            input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") #!!!
            input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
            
          }#If no agree to algorithm's suggestion on annotating to the only subject that is 100#
        }#If only 1 subject is 100%
        
        
        else{#If multiple subjects are 100%
          multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
          colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
          multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus,multiple_sbj$species)),]
          multiple_sbj <- multiple_sbj[,1:2]
          multiple_sbj$comb_name <- rowpaste(multiple_sbj)
          multiple_sbj$pident_max <- NA
          
          #      multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
          #      colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
          if (prod(allduplicates(multiple_sbj$genus))==1) { #If multiple subjects are the same genus
            
            cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
            cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                        magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                        magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
            cat("\n");questions("The algorithm suggests to annotate to genus rank ::");cat("\n")
            cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus)," spp."),30),mean(summary_of_list$pident_max[l])),"::",
                        magic.pident(pr_df(mean(summary_of_list$pident_max[l]),5),mean(summary_of_list$pident_max[l])), 
                        magic.pident(pr_df(mean(summary_of_list$pident[l]),5),mean(summary_of_list$pident[l])))),sep = "\n")
            cat("\n");questions("Do you agree ?")
            slyn <- menu(c("yes","no"))
            if (slyn==1) {#If yes agree to algorithm's suggestion on annotating to the genus rank
              sl_0 <- paste0(paste(unique(multiple_sbj$genus))," spp.",collapse = " | ")
              #         sl <- paste(sl_0," || ", paste(chsl[grep(strsplit(sl_0," ")[[1]][1],chsl)],collapse = " | "))
              sl <- paste(sl_0," || ", paste(multiple_sbj$comb_name,collapse = " | "))
              input$annotated_tax[input$qseqid==i] <- sl
              
              for (ii in length(multiple_sbj$comb_name)) {multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii],summary_of_list$Ssciname))]}
              input$pmatchsel[input$qseqid==i] <- paste(round((multiple_sbj$pident_max/100),3),collapse = " | ")
              
              #          input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[grep(strsplit(sl_0," ")[[1]][1],chsl)]]/100),3),collapse = " | ")
            }#If yes agree to algorithm's suggestion on annotating to the genus rank
            
            else if (slyn==2) { #If no agree to algorithm's suggestion on annotating to the genus rank
              cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
              sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                                       magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                                       magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
                                 "Unidentified","Enter manually"))
              if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")}
              
              #Writing on the fasta file
              input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") 
              input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
              
            }#If no agree to algorithm's suggestion on annotating to the genus rank
          } #If multiple subjects are the same genus
          else if (prod(allduplicates(multiple_sbj$genus))==0) {#If multiple subjects are not the same genus
            cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
            cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                        magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                        magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
            cat("\n");questions("The algorithm suggests to annotate to all subjects with 100% match ::");cat("\n")
            sl <- paste(unique(paste(multiple_sbj$genus,multiple_sbj$species)),collapse = " | ")
            cat(sl);cat("\n")
            cat("\n");questions("Do you agree ?")
            slyn <- menu(c("yes","no"))#If yes agree to algorithm's suggestion on annotating all subjects that are 100
            if (slyn==1) {
              input$annotated_tax[input$qseqid==i] <- sl
              input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%str_trim(unlist(strsplit(sl, "[|]")))]/100),3),collapse = " | ")
            }#If yes agree to algorithm's suggestion on annotating all subjects that are 100
            
            else if (slyn==2) {#If no agree to algorithm's suggestion on annotating all subjects that are 100
              cat("\n");questions("Select the taxa that the query should be annotated or enter the name manually ::");cat("\n")
              sl <- multi.menu(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                                       magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                                       magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident)),
                                 "Unidentified","Enter manually"))
              if(sum(chsl[sl]=="Enter manually")>0){chsl[sl] <- readline("Enter the name of taxa \n(use | between multiple sp) ")}
              
              #Writing on the fasta file
              input$annotated_tax[input$qseqid==i] <- paste(chsl[sl],collapse=" | ") 
              input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%chsl[sl]]/100,3),collapse=" | ") #!!!
            }#If no agree to algorithm's suggestion on annotating all subjects that are 100
          }#If multiple subjects are not the same genus
        }#If multiple subjects are 100%
        
        #Esthetics
        setTxtProgressBar(pb,j);cat("\n") #progress bar
        j <- j+1 #loop mechanism
        write.csv(input,"manually_annotated_function_outcome.csv",row.names = F)
      } #Not skipping method 2
    } #Loop of method 2
  } #Method 2
  
  else if (method==3) { #Method 3
    j=1 #build the loop
    while (j <= length(unique(input$qseqid))){ #Loop of method 3
      pb <- txtProgressBar(1,length(unique(input$qseqid)),style = 3) #progress bar
      i <- unique(input$qseqid)[j] #select the unique query ID
      list_of_query <- input[input$qseqid==i,] #Create a dataframe with all subjects the query blasted towards
      
      #Esthetics of the list
      if (sum(skip==T&!is.na(list_of_query$annotated_tax))>0) {setTxtProgressBar(pb,j);cat("\n");j <- j+1} #Skipping in method 2
      else { #Not skipping method 3
        list_of_query_print <- tibble::tibble(list_of_query);
        list_of_query_print$pident <- magic.pident(list_of_query$pident,list_of_query$pident);
        list_of_query_print$mismatch <- magic.pident(list_of_query$mismatch,list_of_query$pident)
        list_of_query_print$mismatch <- function_col(list_of_query_print$mismatch) 
        list_of_query_print$pident <- function_col(list_of_query_print$pident)
        
        #Print
        if (print_blast==T) {
          cat("\n");questions(paste0("The query \"",i, "\" blast resutls ::"));cat("\n")
          print(list_of_query_print[,c(1:2,5,7:13)]);cat("\n")
        }
        
        #Create a summary list
        summary_of_list <- list_of_query%>% group_by(Ssciname) %>% summarize(pident_max=max(pident), pident=mean(pident)) #Summary of the dataframe
        summary_of_list <- summary_of_list[order(summary_of_list$pident_max,decreasing = T),] #Order the summary based on the %match
        
        chsl <- summary_of_list$Ssciname
        
        match_id <- 100
        l <- summary_of_list$pident_max==match_id
        if (sum(l)==0) {
          l <- summary_of_list$pident_max==max(summary_of_list$pident_max)
          match_id <- max(summary_of_list$pident_max)
        }
        if (sum(l)==1) {#If only 1 subject is 100%
          cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
          cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                      magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                      magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
          cat("\n");questions("The algorithm selected ::")
          cat("\n");cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname[l],30),summary_of_list$pident_max[l]),"::",
                                magic.pident(pr_df(summary_of_list$pident_max[l],5),summary_of_list$pident_max[l]), 
                                magic.pident(pr_df(summary_of_list$pident[l],5),summary_of_list$pident[l]))),sep="\n")
          
          input$annotated_tax[input$qseqid==i] <- paste(summary_of_list$Ssciname[l])
          input$pmatchsel[input$qseqid==i] <- paste(round(summary_of_list$pident_max[summary_of_list$Ssciname%in%summary_of_list$Ssciname[l]]/100,3),collapse=" | ") #Assigns the selected taxa to the input database
          
        }#If only 1 subject is 100%
        
        else{#If multiple subjects are 100%
          multiple_sbj <- unlist_uneven_list(strsplit(str_trim(unlist(strsplit(summary_of_list$Ssciname[l], "[|]"))), "[ ]+"))
          colnames(multiple_sbj) <- c("genus","species",rep("sub.sp",ncol(multiple_sbj)-2))
          multiple_sbj <- multiple_sbj[!duplicated(paste(multiple_sbj$genus,multiple_sbj$species)),]
          multiple_sbj <- multiple_sbj[,1:2]
          multiple_sbj$comb_name <- rowpaste(multiple_sbj)
          multiple_sbj$pident_max <- NA
          if (prod(allduplicates(multiple_sbj$genus))==1) { #If multiple subjects are the same genus
            chsl <- multiple_sbj$comb_name
            cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
            cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,30),summary_of_list$pident_max),"::",
                        magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                        magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
            cat("\n");questions("The algorithm suggests to annotate to genus rank ::");cat("\n")
            cat(c(paste(magic.pident(pr_df(paste0(unique(multiple_sbj$genus)," spp."),30),mean(summary_of_list$pident_max[l])),"::",
                        magic.pident(pr_df(mean(summary_of_list$pident_max[l]),5),mean(summary_of_list$pident_max[l])), 
                        magic.pident(pr_df(mean(summary_of_list$pident[l]),5),mean(summary_of_list$pident[l])))),sep = "\n")
            
            sl_0 <- paste0(paste(unique(multiple_sbj$genus))," spp.",collapse = " | ")
            sl <- paste(sl_0," || ", paste(multiple_sbj$comb_name,collapse = " | "))
            input$annotated_tax[input$qseqid==i] <- sl
            
            for (ii in length(multiple_sbj$comb_name)) {multiple_sbj$pident_max <- summary_of_list$pident_max[min(grep(multiple_sbj$comb_name[ii],summary_of_list$Ssciname))]}
            input$pmatchsel[input$qseqid==i] <- paste(round((multiple_sbj$pident_max/100),3),collapse = " | ")
          } #If multiple subjects are the same genus
          
          else if (prod(allduplicates(multiple_sbj$genus))==0) {#If multiple subjects are not the same genus
            cat("\n");questions("A summary of subjects that your query has been blasted towards ::");cat("\n")
            cat(c(paste(magic.pident(pr_df(summary_of_list$Ssciname,40),summary_of_list$pident_max),"::",
                        magic.pident(pr_df(summary_of_list$pident_max,5),summary_of_list$pident_max), 
                        magic.pident(pr_df(summary_of_list$pident,5),summary_of_list$pident))),sep="\n");cat("\n")
            cat("\n");questions("The algorithm suggests to annotate to all subjects with 100% match ::");cat("\n")
            sl <- paste(unique(paste(multiple_sbj$genus,multiple_sbj$species)),collapse = " | ")
            input$annotated_tax[input$qseqid==i] <- sl
            input$pmatchsel[input$qseqid==i] <- paste(round((summary_of_list$pident_max[summary_of_list$Ssciname%in%str_trim(unlist(strsplit(sl, "[|]")))]/100),3),collapse = " | ")
          }#If multiple subjects are not the same genus
        }#If multiple subjects are 100%
        
        #Esthetics
        setTxtProgressBar(pb,j);cat("\n") #progress bar
        j <- j+1 #loop mechanism
        write.csv(input,"manually_annotated_function_outcome.csv",row.names = F)
      } #Not skipping method 3
    } #Loop of method 3
  } #Method 3
  ret_df <- setNames(data.frame(matrix(NA,length(unique(input$qseqid)),3)),c("qseqid","ssciname","pident"))
  ret_df$qseqid <- c(unique(input$qseqid))
  for (ii in unique(input$qseqid)) {
    ret_df[ret_df$qseqid==ii,2] <- unique(input$annotated_tax[input$qseqid==ii])  
    ret_df[ret_df$qseqid==ii,3] <- unique(input$pmatchsel[input$qseqid==ii])  
  }
  
  return(ret_df)
}#function

