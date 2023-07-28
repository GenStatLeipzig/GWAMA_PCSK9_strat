plotForest3 <- function(mysnpstats,data) {
  # debug
  # mysnpstats = ToDoList[1,]
  # data = copy(myTab)

  mysnpstatsm = copy(data)
  mysnpstatsm = mysnpstatsm[markername == mysnpstats$SNP,]
  mysnpstatsm[,table(phenotype)]
  setorder(mysnpstatsm,betaFEM) 
  
  # calc ranges to plot, default is not nice, since reference line drops out of sight for large values
  mysnpstatsm[,lowerCI95 := betaFEM-1.96*seFEM]
  mysnpstatsm[,upperCI95 := betaFEM+1.96*seFEM]
  xmin = min(c(0,mysnpstatsm$lowerCI95, mysnpstatsm$upperCI95), na.rm = T)
  xmax = max(c(0,mysnpstatsm$lowerCI95, mysnpstatsm$upperCI95), na.rm = T)
  xrange = xmax - xmin
  xmin_margin = xmin - 1.5*xrange
  xmax_margin = xmax + 1.5*xrange
  
  # start plot
  par(mar=c(4,6,0,4))
  par(font=1)
  
  myPheno = mysnpstatsm$phenotype
  myGene = mysnpstats$candidateGene
  mySNP = gsub(":.*","",mysnpstats$SNP)
  meanEAF = round(mysnpstatsm[,mean(nWeightedEAF)],3)
  meanInfo = round(mysnpstatsm[,mean(nWeightedInfoScore)],3)
  
  dets = forest(x=mysnpstatsm$betaFEM, 
                sei = mysnpstatsm$seFEM,
                #xlab=paste0(mySNP," (EAF=",meanEAF,"; info=",meanInfo,"; gene=",myGene,")"),
                xlab = "",
                showweights=T, 
                slab = mysnpstatsm$phenotype,
                ylim = c(0, nrow(mysnpstatsm)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin, xmax))
  myText = paste0(mySNP," (EAF=",meanEAF,"; info=",meanInfo,"; gene=",myGene,")")
  par(font=4)
  text(min(dets$xlim), max(dets$ylim-0.5),  myText , pos=4)
  par(font=2)
  text(min(dets$xlim), max(dets$ylim-1.5), "PCSK9 subgroups" , pos=4)
  text(max(dets$xlim), max(dets$ylim-1.5), "Weight         Effect (95% CI)", pos=2)
  
}
