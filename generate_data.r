library(cgdsr)
library(survival)
library(survminer)
source("/home/gabriela/Documents/code/Cox/GSA.read.gmt.r")
source("/home/gabriela/Documents/code/Cox/pathway_analysis.r")

mycgds = CGDS("http://www.cbioportal.org/")

for(ii in mySelectedStudies){ #get this by running gettingTCGAsamplesWithExpressio...
gflag = TRUE #FALSE  # TRUE if genes should be grouped into pathways
alp = 1  # elastic net
# for(K in c(10,20,30,40,50)){
K <- 1   # nr of cross-validation sets
nmax = 1 #30   # number of repetitions
adapt = FALSE  # adaptive choice of low/high risk split
paths <- NULL
path.names <- NULL

aas = read.csv("/home/gabriela/Documents/code/Cox/pathways/pathways_no_dupl.csv",header = TRUE)
aas = aas[,-1]
pathway_names = colnames(aas)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudies = getCancerStudies(mycgds)[,1]
mycancerstudy=getCancerStudies(mycgds)[ii,1]
studyName=getCancerStudies(mycgds)[ii,2]

## Are there gene expression data?
caselist = getCaseLists(mycgds,mycancerstudy)[,2];	# all case expressions
tt=which(caselist=="Tumor Samples with mRNA data (RNA Seq V2)" | caselist=="Samples with mRNA data (RNA Seq V2)" )
mycaselist = getCaseLists(mycgds,mycancerstudy)[tt,1]  #[1,1]
# expressionPerTumor=getProfileData(mycgds,myListOfHousekeeping,mycaselist,mycaselist)

#############################################################################################
###### Clinical data ########################################################################
#############################################################################################
clinicalPerTumor=getClinicalData(mycgds,mycaselist)

# extract survival times and status
os_months_index = which(names(clinicalPerTumor)=="OS_MONTHS")  # find the index of time column
os_status_index = which(names(clinicalPerTumor)=="OS_STATUS")  # find the index of status column

# select clinical covariates
covariates <- NULL #c(1,2,3) #c(1,3,4,8,21,40)  # choose covariates to include in the study
x_first = data.frame(clinicalPerTumor[,c(os_months_index, os_status_index)]) # time & status
if(!is.null(covariates)){  # if some clinical covariates were chosen
  x_long = data.frame(clinicalPerTumor[,covariates])

  ## all categorical into one-hot
  dmy <- dummyVars(" ~ .", data = x_long)   # dummify the data
  aas <- data.frame(predict(dmy, newdata = x_long))
  names <- colnames(aas)  # save the names
  xx1 <- cbind(x_first,aas)  # combine
  # xx <- model.matrix( ~ ., xx)[,-1]  # dataframe -> matrix
  ne_clinical <- ncol(xx1)-2  # nr of clinical variables (w/o status and time)

  # define the groups
  group <- 1:ne_clinical
} else{ # no clinical covariates
  x_long = NULL
  names = NULL
  xx1 = x_first
  group = NULL
  ne_clinical <- 0
}
#################################################################################################
############### Gene expression data ############################################################
#################################################################################################

# append expression data
xx2 = NULL  # initialize
# for(i in 1:length(paths)){
for(i in 1:length(aas)){
  # expr = getProfileData(mycgds,paths[[i]],mycaselist,mycaselist)  # Expression per tumor per path
  expr = getProfileData(mycgds,aas[,i],mycaselist,mycaselist)
  ne_genes = dim(expr)[2]  # how many genes per path
  col2remove <- NULL
  # remove genes with NA data (if the whole column is NA)
  for(jj in 1:ne_genes){
        if(any(is.na(expr[,jj])) || all(expr[,jj]==0)){
          col2remove <- c(col2remove,jj) # mark columns for removal
        }
      }
  if(!is.null(col2remove)){  # if there were any NaN or 0 columns found
    expr <- expr[-col2remove]  # delete NaN and 0 columns
  }
  ne_genes_new <- ne_genes-length(col2remove)  # how many genes left after removing NaN columns
    if(i==1){ # in the first pathway, xx2 is not appendable
      xx2 = expr
    } else if(length(expr)==0){} # don't try to append if expr is empty
    else {
      xx2 = data.frame(cbind(xx2,expr)) # append expression data to design matrix
    }
  # xx = data.frame(c(xx,expr)) # append expression data to design matrix
  group = c(group,rep.int(ne_clinical+i,ne_genes_new))  # append the group number
}

# append epression data # BUT SORT THEM FIRST
xx1_sorted<-xx1[ order(rownames(xx1)),]  # clinical
xx2_sorted<-xx2[ order(rownames(xx2)),]  # gene expression
# check if the patient set is the same
stopifnot(all(rownames(xx1_sorted)==rownames(xx2_sorted)))
xx_sorted = data.frame(cbind(xx1_sorted,xx2_sorted))


names <- c(names,path.names)  # add names to the list
if(!gflag){
  group = 1:(ncol(xx_sorted)-2)  # rewrite if the genes are not grouped
}
#################################################################################################
############### Comb and fit the model ##########################################################
#################################################################################################

## remove NA and other missing data
ne_long <- dim(xx_sorted)[1]  # how many samples before pruning
ind_to_remove2 = which(xx_sorted[,1] == 0)   # find cases where t = 0
if(!length(ind_to_remove2)==0){ # no elements to remove
  xx_sorted = xx_sorted[-ind_to_remove2,]   # remove t=0 elements from design matrix
}
xx_sorted = na.omit(xx_sorted)  # remove NA elements
ne_short <- dim(xx_sorted)[1] # nr of samples after pruning
warning(ne_long-ne_short," NA elements removed")

# split the big matrix into x and y data
os_months = xx_sorted[,1]
os_status = (xx_sorted[,2]=="1:DECEASED")  # find indices of all deceased # previously just DECEASED
x = xx_sorted[,3:ncol(xx_sorted)]

# create the object
surv_object <- Surv(time = os_months, event = os_status)
  
xxx <- model.matrix( ~ ., x)[,-1]     		# transform x into a matrix

ctype = sub("\\_.*", "", mycancerstudy) # study abbreviation
print(paste("This",studyName,"is going to be legen...", toString(jj),"/",toString(nmax)))
write.csv(group,paste0('/home/gabriela/Documents/code/Cox/surv_files/',ctype,'_groups.csv'))
write.csv(os_months,paste0('/home/gabriela/Documents/code/Cox/surv_files/',ctype,'_os_months.csv'))
write.csv(os_status,paste0('/home/gabriela/Documents/code/Cox/surv_files/',ctype,'_os_status.csv'))
write.csv(xxx, paste0('/home/gabriela/Documents/code/Cox/surv_files/',ctype,'_x.csv'))
}
