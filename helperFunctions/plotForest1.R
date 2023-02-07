plotForest1 <- function(mysnpstats,data) {
  # debug
  # mysnpstats = tab1[1,]
  # data = copy(tab2)

  mysnpstatsm = copy(data)
  mysnpstatsm = mysnpstatsm[markername == mysnpstats$markername,]
  mysnpstatsm = mysnpstatsm[pheno == mysnpstats$pheno,]
  
  mysnpstatsm2 = dcast.data.table(mysnpstatsm, markername +pheno    +betaFEM   +   seFEM   +     pFEM+ cochQ + cochQpval + I2   + betaREM     + seREM +     pREM + cohort~stats, value.var = "value")
  mysnpstatsm2[,table(cohort)]
  
  # calc ranges to plot, default is not nice, since reference line drops out of sight for large values
  mysnpstatsm2[,lowerCI95 := beta-1.96*se]
  mysnpstatsm2[,upperCI95 := beta+1.96*se]
  xmin = min(c(0,mysnpstatsm2$lowerCI95, mysnpstatsm2$upperCI95), na.rm = T)
  xmax = max(c(0,mysnpstatsm2$lowerCI95, mysnpstatsm2$upperCI95), na.rm = T)
  xrange = xmax - xmin
  xmin_margin = xmin - 1.5*xrange
  xmax_margin = xmax + 1.5*xrange
  
  # start plot
  par(mar=c(5,6,0,4))
  par(font=1)
  
  myPheno = gsub("PCSK9_","",mysnpstats$pheno)
  myPheno = gsub("_"," ",myPheno)
  
  dets = forest(x=mysnpstatsm2$beta, 
                sei = mysnpstatsm2$se,
                xlab=myPheno,
                showweights=T, 
                slab = mysnpstatsm2$cohort,
                ylim = c(-2, nrow(mysnpstatsm2)+3) , 
                cex =1, 
                xlim = c(xmin_margin, xmax_margin), 
                alim = c(xmin, xmax))
  par(font=4)
  text(min(dets$xlim), max(dets$ylim-1.5),  bquote(paste("Study (",I^2,"=",.(round(mysnpstats$I2*100)), "%)")) , pos=4)
  text(max(dets$xlim), max(dets$ylim-1.5), "Weight         Effect (95% CI)", pos=2)
  par(font=2)
  addpoly(x = c(mysnpstats$betaFEM, mysnpstats$betaREM),
          sei =  c(mysnpstats$seFEM, mysnpstats$seREM),
          row=c(-0.5,-1.5),
          cex=1,
          mlab=c(paste0("Fixed Effect (p=",signif(mysnpstats$pFEM, 2), ")") ,
                 paste0("Random Effect (p=",signif(mysnpstats$pREM, 2), ")")))

}
