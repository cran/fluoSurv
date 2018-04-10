## these routines are intended to read fluorescence data as produced by a Biotek plate reader ##

## reads a bloc of data
read.block.base <- function(txt,puits=NULL) {
  pos.ligne <- 2:length(txt)
  pos.ligne <- pos.ligne[grep("^[A-H]\t",txt[pos.ligne])]
  if(is.null(puits)) {
                                        #reads column number
    col <- strsplit(txt[1],"\t")[[1]]
    col <- col[-1]
                                        #reads line number
    ligne <- substr(txt[pos.ligne],1,1)

    puits   <- paste(ligne[sort(rep(1:length(ligne),length(col)))],rep(col,length(ligne)),sep="")
  }
                                        #reads measurements
  zz <- strsplit(gsub(",",".",txt[pos.ligne]),"\t")
  conversion <- function(y) ifelse(y=="OVRFLW",NA,y)
  zz <- sapply(zz,function(x) conversion(x[-c(1,length(x))]))
                                        # first element is line number, last element is read ID
                                        # the two elements must be discarded
  return(data.frame(well=puits,value=as.numeric(zz)))
}

## reads bloc of data and ID
read.block <- function(pos1,pos2,l,puits=NULL) {
  #identifiant
  ID_read <- gsub(",","_",gsub("^Read ([0-9]*):(.*)$","\\1_\\2",l[pos1]))
  #lecture des mesures
  dtmp <- read.block.base(l[(pos1+1):pos2],puits)
  dtmp$ID_read <- ID_read
  return(dtmp)
}

## reads a file with a single measurement
read.spectro <- function(nom,path) {
  cat("Chargement du fichier",nom,"\n")
  l <- readLines(paste(path,nom,sep="/"))

  ## position des lectures dans le fichier
  pos <- grep("^Read [0-9]*:.*",l)
  if(length(pos)<=0) return(NULL)
  next.pos <- c(pos[-1]-1,length(l))

  d <- read.block(pos[1],next.pos[1],l)
  liste.puits <- as.character(d$puits)
  taille.bloc <- nrow(d)
  pos.bloc <- taille.bloc
  d <- rbind(d,data.frame(well=rep(liste.puits,length(pos)-1),value=NA,ID_read=NA))
  for(i in 2:length(pos)) {
    d[pos.bloc+1:taille.bloc,] <- read.block(pos[i],next.pos[i],l,liste.puits)
    pos.bloc <- pos.bloc+taille.bloc
  }

  nom.sortie <- gsub("^(.*)txt$","table_\\1csv",nom)
  cat("Save data in",nom.sortie,"\n")
  write.table(file=nom.sortie,d)
}

## reads time values
strtotime <- function(str) {
  str <- gsub("^.*\\((.*)\\)$","\\1",str)
  v <- strsplit(str,":")[[1]]
  sum(as.numeric(v)/c(1,60,3600))
}

## reads a bloc in a kinetic file
read.block.kinetic <- function(pos1,pos2,l,liste.puits=NULL,readTime=T) {
  t <- NA
  if(readTime) {
    #lecture du temps
    t <- strtotime(as.vector(l)[pos1]);
  }
  #reads number
  num <- as.numeric(gsub("^Kinetic read ([0-9]*) \\(.*$","\\1",l[pos1]))
  # reads ID
  ID_read <- gsub(",","_",gsub("^Read ([0-9]*):(.*)$","\\1_\\2",l[pos1-1]))
  #reads gain and wavelengths
  read <- strsplit(ID_read,"_")[[1]]
  #reads measurement
  dtmp <- read.block.base(l[(pos1+1):pos2],puits=liste.puits)

  dtmp$t <- t
  dtmp$num <- num
  dtmp$read <- ifelse(length(read)>=1,as.numeric(read[1]),NA)
  dtmp$exc  <- ifelse(length(read)>=2,as.numeric(read[2]),NA)
  dtmp$em   <- ifelse(length(read)>=3,as.numeric(gsub("[^0-9].*$","",read[3])),NA) # gsub serves when different gain values are used for the same pairs of
                                                                                   # wavelengths. Wevelength is then followed by a number in [!
  dtmp$ID_read <- ID_read
  return(dtmp)
}

## reads a whole kinetic file

#' Reads a kinetic file, as produced by a Biotek plate reader.
#'
#' @param name The name of the file to be read
#' @param path The path where the file is to be found
#' @param readTime Should time data be read?
#' @param saveData Should the resulting \code{data.frame} be saved?
#'
#' @return Returns a \code{data.frame} if \code{saveData} is set to \code{FALSE}.
#' If \code{saveData} is set to \code{TRUE}, the \code{data.frame} is saved and the
#' file name is returned.
#'
#' @export
#'
#' @examples
#' ## reads data. Warning: files are large, and this operation takes time!
#' d <- read.kinetic("kinetics_xenorhabdus_galleria.txt",
#'                      path=system.file('extdata', package = 'fluoSurv'),
#'                      saveData=FALSE)
#' str(d)
#'
#' ## saveData should rather be set to TRUE so that converted data are saved
#' ## in a csv file and can be re-used later on.
#'
read.kinetic <- function(name,path=NULL,readTime=TRUE,saveData=TRUE) {
  cat("Loads kinetic file",name,"\n")
  outfile <- gsub("^(.*)txt$","table_\\1csv",name)
  if(!is.null(path)) name <- paste(path,name,sep="/")
  l <- readLines(name)

  ## position des lectures dans le fichier
  pos <- grep("^Kinetic read",l)
  cat(length(pos)," blocks detected in ",name,"\n")

  next.pos <- c(pos[-1]-1,length(l))

  cat(" [0] .")
  d <- read.block.kinetic(pos[1],next.pos[1],l,readTime=readTime)
  liste.puits <- as.character(d$well)
  taille.bloc <- nrow(d)
  pos.bloc <- taille.bloc
  dtmp <- data.frame(well=rep(liste.puits,length(pos)),value=NA,t=NA,num=NA,read=NA,exc=NA,em=NA,ID_read=NA)
  d <- rbind(d,dtmp)
  for(i in 2:length(pos)) {
    cat(".")
    if(i - 50*(i %/% 50) == 0) cat("\n [",i,"]")
    d[pos.bloc+1:taille.bloc,] <- read.block.kinetic(pos[i],next.pos[i]-1,l,liste.puits,readTime)
    pos.bloc <- pos.bloc+taille.bloc
  }

  if(!saveData) return(d)


  cat("Save tabular data in",outfile,"\n")
  write.table(file=outfile,d,sep=";")
  return(outfile)
}


