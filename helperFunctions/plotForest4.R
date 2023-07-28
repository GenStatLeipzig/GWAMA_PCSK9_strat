plotForest4 <- function(mysnpstats,data) {
  # debug
  # mysnpstats = ToDoList[1,]
  # data = copy(myTab[chr==1])

  mysnpstatsm = copy(data)
  #mysnpstatsm = mysnpstatsm[markername == mysnpstats$SNP,]
  mysnpstatsm[,table(phenotype)]
  setorder(mysnpstatsm,pos,betaFEM) 
  
  # calc ranges to plot, default is not nice, since reference line drops out of sight for large values
  mysnpstatsm[,lowerCI95 := betaFEM-1.96*seFEM]
  mysnpstatsm[,upperCI95 := betaFEM+1.96*seFEM]
  xmin = min(c(0,mysnpstatsm$lowerCI95, mysnpstatsm$upperCI95), na.rm = T)
  xmax = max(c(0,mysnpstatsm$lowerCI95, mysnpstatsm$upperCI95), na.rm = T)
  xrange = xmax - xmin
  xmin_margin = xmin - 0.5*xrange
  xmax_margin = xmax + 0.5*xrange
  
  # start plot
  par(mar=c(5,6,0,4))
  par(font=1)
  
  myGene = "PCSK9"
  mySNP = unique(gsub(":.*","",mysnpstatsm$markername))
  meanEAF = mysnpstatsm[,mean(nWeightedEAF),by="markername"]
  setnames(meanEAF,"V1","EAF")
  meanEAF[,EAF := round(EAF,3)]
  meanInfo = mysnpstatsm[,mean(nWeightedInfoScore),by="markername"]
  setnames(meanInfo,"V1","info")
  meanInfo[,info := round(info,3)]
  dummy = copy(mysnpstatsm)
  dummy = dummy[c(1,9,17),]
  dummy[,betaFEM :=NA]
  dummy[,seFEM :=NA]
  dummy[,phenotype := paste0(gsub(":.*","",meanEAF$markername), " (EAF=",meanEAF$EAF,"; info=",meanInfo$info,")")]
  dummy[,pos := pos - 1]
  dummy2 = copy(dummy)
  dummy2[,pos := pos - 1]
  dummy2 = dummy2[-1,]
  dummy2[,phenotype := NA]
  
  mysnpstatsm = rbind(mysnpstatsm,dummy,dummy2)
  setorder(mysnpstatsm,pos, betaFEM) 
  
  
  options(na.action = "na.pass")
  dets = forest(x=mysnpstatsm$betaFEM, 
                sei = mysnpstatsm$seFEM,
                xlab="independent SNP effects at PCSK9 gene locus",
                showweights=F, 
                slab = mysnpstatsm$phenotype,
                ylim = c(0, nrow(mysnpstatsm)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin, xmax))
  par(font=4)
  text(min(dets$xlim), max(dets$ylim-1.5),  bquote(paste("PCSK9 Subgroups")) , pos=4)
  text(max(dets$xlim), max(dets$ylim-1.5), "         Effect (95% CI)", pos=2)
  options(na.action = "na.omit")
  
}
