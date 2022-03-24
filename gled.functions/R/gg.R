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
#' @return otu table with only presence / absence
#' @export
#' @examples
#' otu<-pa_conversion(otu)


pa_conversion<-function(otu){
  otu[otu>0]<-1
  return(as.data.frame(otu))
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
require(plyr)
require(stringr)
require(ggplot2)
motudiv<-function(taxa, rank, logic=1){
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

