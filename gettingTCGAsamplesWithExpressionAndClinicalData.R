# install.packages("cgdsr")

library(cgdsr)

# outputDir="C:/WorkFolder/TCGA/"
outputDir="/home/gabriela/Documents/code/Cox/"
setwd(outputDir)

# Create CGDS object
# mycgds = CGDS("http://www.cbioportal.org/public-portal/")
mycgds = CGDS("http://www.cbioportal.org/")
getCancerStudies(mycgds)
# Get available case lists (collection of samples) for a given cancer study
mycancerstudies = getCancerStudies(mycgds)[,1]

#create baseline for expression:
myListOfHousekeeping=c("GAPDH","ACTB","GUSB","B2M","HMBS","HPRT1","RPL13A","SDHA","TBP","UBC","YWHAZ")

samplesWithData=0
myBaseline=NULL
sampleNames=NULL
mySelectedStudies=NULL
NumberOfSamples=NULL
hasToLog=NULL #whether we need to transform in log
typeOfExprProtocol=NULL #RNA-seq V2, RNA-seq or microarrays
for (i in c(1:length(mycancerstudies))) { #239 studies
  mycancerstudy=getCancerStudies(mycgds)[i,1]
  studyName=getCancerStudies(mycgds)[i,2]
  print(paste0("Sample    ",studyName))
  
  tt=which(getCaseLists(mycgds,mycancerstudy)[,2]=="Tumor Samples with mRNA data (RNA Seq V2)" | getCaseLists(mycgds,mycancerstudy)[,2]=="Samples with mRNA data (RNA Seq V2)" )
  typeOfExpressionThisCancer="RNA-seq V2"
  # if (length(tt)==0) {
  #   tt=which(getCaseLists(mycgds,mycancerstudy)[,2]=="Tumor Samples with mRNA data (RNA Seq)" | getCaseLists(mycgds,mycancerstudy)[,2]=="Samples with mRNA data (RNA Seq)")
  #   typeOfExpressionThisCancer="RNA-seq"
  # }
  
  # if (length(tt)==0) {
  #   tt=which(getCaseLists(mycgds,mycancerstudy)[,2]=="Tumors with mRNA data (Agilent microarray)" | getCaseLists(mycgds,mycancerstudy)[,2]=="Samples with mRNA data (Agilent microarray)")
  #   typeOfExpressionThisCancer="Microarray"
  # }
  
  if (length(tt)) {
    
    print("Found expression data")
    mycaselist = getCaseLists(mycgds,mycancerstudy)[tt,1]
    
    result=tryCatch({
      clinicalPerTumor=getClinicalData(mycgds,mycaselist)
      1
    }, error = function(e) {
      0
    })
    
    if (result==1) {
      print("Found clinical data")
      if ("OS_STATUS" %in% names(clinicalPerTumor) && "OS_MONTHS" %in% names(clinicalPerTumor) ) {
        print("Found survival data")
        
        print(paste("Will use for expression: ",getCaseLists(mycgds,mycancerstudy)[tt,2]))
        expressionPerTumor=getProfileData(mycgds,myListOfHousekeeping,mycaselist,mycaselist) 
        if (!is.na(mean(colMeans(expressionPerTumor,na.rm = T),na.rm = T))) {
          
          myMin=min(expressionPerTumor,na.rm = T)
          myMax=max(expressionPerTumor,na.rm = T)
          
          typeOfExprProtocol=c(typeOfExprProtocol,typeOfExpressionThisCancer)
          
          if (!is.na(myMin)) {
            if (myMin>=0 && myMax>30) {
              expressionPerTumor=log(expressionPerTumor+1,base = 2)
              print ("Transforming Expression Values to LOG2")
              hasToLog=c(hasToLog,TRUE)
            }else {hasToLog=c(hasToLog,FALSE)}
            print(paste0('Success! - found all information for ',nrow(expressionPerTumor)," tumors"))
            myBaseline=c(myBaseline,mean(colMeans(expressionPerTumor,na.rm = T),na.rm = T))
            mySelectedStudies=c(mySelectedStudies,i)
            sampleNames=c(sampleNames,studyName)
            samplesWithData=samplesWithData+1
            NumberOfSamples=c(NumberOfSamples,nrow(expressionPerTumor))
          }
        }
      }
    }
  }
  
}

View(cbind(sampleNames,NumberOfSamples,myBaseline,hasToLog,typeOfExprProtocol))

#remove smaller studies for the same cancer:
studyNamesShort=sub(" [(].*", "", sampleNames, perl = T)
indToRemove=tapply(c(1:length(studyNamesShort)), studyNamesShort, function(x){
  if(length(x)>1) {
    # print(x) 
    # maxToKeep=which.max(NumberOfSamples[x])
    maxToKeep=which(NumberOfSamples[x] == max(NumberOfSamples[x], na.rm = TRUE))[length(which(NumberOfSamples[x] == max(NumberOfSamples[x], na.rm = TRUE)))]
    x[-maxToKeep]
  }
})
indToRemove=unlist(indToRemove)

# sum(!is.na(myBaseline[-indToRemove]))
# sum(!is.na(myBaseline))

mySelectedStudies=mySelectedStudies[-indToRemove]
samplesWithData=length(mySelectedStudies)
myBaseline=myBaseline[-indToRemove]
sampleNames=sampleNames[-indToRemove]
NumberOfSamples=NumberOfSamples[-indToRemove]
hasToLog=hasToLog[-indToRemove]
typeOfExprProtocol=typeOfExprProtocol[-indToRemove]
sum(NumberOfSamples)

View(cbind(sampleNames,NumberOfSamples,myBaseline,hasToLog,typeOfExprProtocol))


