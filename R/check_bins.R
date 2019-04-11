check_bedtools = function(path = NULL){

  if(is.null(path)){
    check = as.character(Sys.which(names = 'bedtools'))[1]
  }else{
    check = paste0(path, "/bedtools")
    if(!file.exists(check)){
      stop("Could not locate bedtools at " , check ,"\nDownload bedTools: <http://code.google.com/p/bedtools/>")
    }
  }

  if(check != ""){
    return(check)
  }else{
    stop("Could not locate bedtools at " , check ,"\nDownload bedTools: <http://code.google.com/p/bedtools/>")
  }
}

#---------------------------------------------------------------------------------------------------------------

check_bedGraphToBigWig = function(path = NULL){

  if(is.null(path)){
    check = as.character(Sys.which(names = 'bedGraphToBigWig'))[1]
  }else{
    check = paste0(path, "/bedGraphToBigWig")
    if(!file.exists(check)){
      stop("Could not locate bedGraphToBigWig at " , check ,"\nDownload: <http://hgdownload.cse.ucsc.edu/admin/exe/>")
    }
  }

  if(check != ""){
    return(check)
  }else{
    stop("Could not locate bedGraphToBigWig at " , check ,"\nDownload: <http://hgdownload.cse.ucsc.edu/admin/exe/>")
  }
}


#---------------------------------------------------------------------------------------------------------------

check_bedClip = function(path = NULL){

  if(is.null(path)){
    check = as.character(Sys.which(names = 'bedClip'))[1]
  }else{
    check = paste0(path, "/bedClip")
    if(!file.exists(check)){
      stop("Could not locate bedClip at " , check ,"\nDownload: <http://hgdownload.cse.ucsc.edu/admin/exe/>")
    }
  }

  if(check != ""){
    return(check)
  }else{
    stop("Could not locate bedClip at " , check ,"\nDownload: <http://hgdownload.cse.ucsc.edu/admin/exe/>")
  }
}
