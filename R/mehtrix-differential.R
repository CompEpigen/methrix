#' Differential Methylation Calling with DSS
#' @details Wrapper function for differential methylation calling using \code{\link{DSS}}
#' @param m Input \code{\link{methrix}} object.
#' @param pheno Column name of colData(m). The values of this column will be used as classes to call differential methylation.
#' @param smooth Should the methylation values be smoothed?
#' @param s_span Smoothing span. See \code{\link{DMLtest}}. Default 500.
#' @param delta Methylation difference to call a DML/DMR significant. Default 0.1.
#' @param p P-value cut-off. Default 0.05.
#' @param minlen Minimal length of a DMR. Default 50 bp.
#' @param minCG Minimal number of CpGs in a DMR. Default 3.
#' @param mergeDis Maximal distance of two consecutive DMRs to be merged. Default 500 bp.
#' @param pct.sig Minimum percentage of significant CpGs in a DMR. Default 0.5.
#' @param parallel Execute job parallel? Just possible if chrWise=T. Default FALSE.
#' @param chrWise Compute job chromosome-wise and compute FDR adjusted p-values afterwards? Default FALSE.
#' @param ncore Number of cores if executed parallel. Default 1.
#'
#' @import DSS
#' @import bsseq
#' @import doParallel
#' @import foreach
#' @return
#' @export
#'
#' @examples
methrix_differential <- function(m,pheno=NULL,smooth=FALSE,s_span=500,delta=0.1,p=0.05,
                                    minlen=50,minCG=3,mergeDis=500,pct.sig=0.5,parallel=F,chrWise=F,ncore=1){

  ## check if the classes are available in the selected phenotype column
  if(chrWise==F&parallel==T){stop("Parallel computing can just be executed in a chromosome wise manner.")}
  if(is.null(pheno)){"You forgot to specify a phenotype column for DMR calling."}
  if(sum(pheno%in%colnames(colData(m)))==0){"You forgot to specify a phenotype column for DMR calling."}
  if(length(levels(as.factor(colData(m)[pheno][,1])))>2){"Methrix currently only supports DMR calling with DSS for two classes. Please choose a phenotype column with just two classes."}
  work_fac <- as.factor(colData(m)[pheno][,1])

  ## define the classes
  class1 <- rownames(colData(m)[which(work_fac==levels(work_fac)[1]),])
  class2 <- rownames(colData(m)[which(work_fac==levels(work_fac)[2]),])
  cat(paste0("The following sample is included in class I: ",class1,"\n"))
  cat(paste0("The following sample is included in class II: ",class2,"\n"))


  ## calculate the DMLs not in a chr wise fashion
  if(chrWise==F&parallel==F){
    ## change the format of the methrix object to a bsseq object
    methrix_bs <-  methrix::methrix2bsseq(m = m)
    dml_test <- DSS::DMLtest(methrix_bs, group1=class1,group2=class2,smoothing=smooth,smoothing.span=s_span)
   # print(head(dml_test))
  }else if(parallel==T){
    chroms <- levels(as.factor(m@elementMetadata@listData$chr))

    #Set-up the cluster
    doParallel::registerDoParallel(cores=ncore)
    cat(paste0("You are now working on ",foreach::getDoParWorkers()," cores.","\n"))
    dml_classI_classII <- foreach::foreach(i=chroms,
                                           .final = function(x) setNames(x, chroms))%dopar%{
      cat(paste0("Currently processing: ",i,"\n"))
      methrix_work <- methrix::subset_methrix(m,contigs = i)
      methrix_bs <-  methrix::methrix2bsseq(m = methrix_work)
      DSS::DMLtest(methrix_bs,group1 = class1, group2 = class2, smoothing=smooth,smoothing.span=s_span)
                                           }
    # combine the DMLs and adjust the p-values
    dml_test<-do.call("rbind",dml_classI_classII)
    dml_test$fdr<-p.adjust(dml_test$pval,method="BH")
  }else if(chrWise==T&parallel==F){
    chroms <- levels(as.factor(m@elementMetadata@listData$chr))
    dml_classI_classII <- list()
    for(i in chroms){
     cat(paste0("Currently processing: ",i,"\n"))
     methrix_work <- methrix::subset_methrix(m,contigs = i)
     methrix_bs <-  methrix::methrix2bsseq(m = methrix_work)
     dml_classI_classII[[i]] <- DSS::DMLtest(methrix_bs,group1 = class1, group2  = class2,smoothing=smooth,smoothing.span=s_span)
                                         }
    # combine the DMLs and adjust the p-values
    dml_test<-do.call("rbind",dml_classI_classII)
    dml_test$fdr<-p.adjust(dml_test$pval,method="BH")
    }

    # call significant DMLs
    sigDML <- DSS::callDML(dml_test, delta=delta, p.threshold=p)
    #print(head(sigDML))
    # call DMRs
    DMRs <- DSS::callDMR(dml_test,delta=delta,p.threshold = p,minlen = minlen,minCG = minCG,dis.merge = mergeDis,pct.sig = pct.sig)

    # return a list with DMLs and DMRs
    return(list("DMLs"=dml_test,"Significant_DMLs"=sigDML,"DMRs"=DMRs))
}






