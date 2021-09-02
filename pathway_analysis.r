analyze_pathways <- function(b,path.names,mycancerstudy,studyName){
library(RColorBrewer)  # color palettes
  # brewer color scales
colscales = c("Blues","BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd",
              "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")
sc_choice = sample(1:length(colscales),1) # make a random color choice

b_length = length(b)
b_unique_unsorted = NULL
for(ii in 1:b_length){
  b_unique_unsorted = c(b_unique_unsorted,b[[ii]])
}
b_unique = sort(unique(b_unique_unsorted))   # all unique pathways
nb = length(b_unique)  # nr of unique elements
b_abs = rep(0,nb)   # absolute appearance of pathway
for(ii in 1:length(b_unique_unsorted)){ # count occurances
  ind = which(b_unique == b_unique_unsorted[ii])
  b_abs[ind] = b_abs[ind] + 1
}
pnames = path.names[b_unique]  # extract pathway names

# appearance frequency
b_mat = matrix(0L,nb,nb)  # correlation matrix
# b_abs = rep(0,nb)   # absolute appearance of pathway
for(ii in 1:b_length){
  pathways = unique(b[[ii]])   # one pathways set at a time # there shouldnt be repetititons but just in case
  nnp = length(pathways)   # nr of pathways in CV set
  if(!nnp==1){  # we need at least two pathways
    for(jj in 1:(length(pathways)-1)){
      ind1 = which(pathways[jj]==b_unique)   # index of the pathway in the matrix
      # b_abs[ind1] = b_abs[ind1] + 1   # raise the absolute count
      for(kk in (jj+1):length(pathways)){  # update the correlation with +1 if i and j are in the same set
        ind2 = which(pathways[kk]==b_unique)   # index of the pathway in the matrix
        b_mat[min(ind1,ind2),max(ind1,ind2)] = b_mat[min(ind1,ind2),max(ind1,ind2)] + 1 # saved consistently in the same triangular matrix
      }
    }
  }
}
# rescale the matrix with total appearance number
b_mat_rescale = matrix(0L,nb,nb)
rownames(b_mat_rescale) = pnames
colnames(b_mat_rescale) = pnames
for(ii in 1:nb){
  for(jj in 1:nb){
    b_mat_rescale[ii,jj] = (b_mat[ii,jj])^2/(b_abs[jj]*b_abs[ii])  # normalize
  }
}

# heatmap
cc <- colorRampPalette(brewer.pal(nb, colscales[sc_choice]))(nb) # interpolate in case not enough colors in scale
s1 = paste0(mycancerstudy,"_K",toString(K),"_heatmap.pdf")
pdf(s1)
heatmap(b_mat_rescale + t(b_mat_rescale),Rowv = NA,Colv = NA, scale = "none",
        main = "Pathway correlation", xlab = studyName, margins = c(15,15),
        col = cc)
dev.off()
# image(t(b_mat_rescale)[ncol(b_mat_rescale):1,])

# Barplot
b_abs_mat = matrix(0L,nb,1)  # relative occurance matrix
b_abs_mat[1:nb,] = 100*b_abs/b_length
rownames(b_abs_mat) = pnames
s2 = paste0(mycancerstudy,"_K",toString(K),"_barplot.pdf")
pdf(s2)
par(mai=c(1,2.5,1,1))
barplot(t(b_abs_mat),#xlab = "Pathways", y
        xlab = "Relative occurrence (%)",#cex.names=.5,
        main="Relative occurrence of pathways", horiz=TRUE, las = 1, # labels and bars running vertically
        col = cc[ceiling(nb/2)])#,col = "#56B4E9") # plot in average heatmap color
dev.off()
}
