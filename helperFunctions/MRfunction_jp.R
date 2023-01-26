MRfunction_jp<-function(betaX,seX,betaY,seY){
  betaIV<-betaY/betaX
  se.ratio.st<- seY/sqrt(betaX^2)
  se.ratio.nd<- sqrt(seY^2/betaX^2 + betaY^2*seX^2/betaX^4)
  p1<-2*pnorm(-abs(betaIV/se.ratio.st))
  p2<-2*pnorm(-abs(betaIV/se.ratio.nd))
  res<-data.table(
    beta_IV=betaIV,
    se_IV1=se.ratio.st,
    se_IV2=se.ratio.nd,
    p_IV1=p1,
    p_IV2=p2)
  return(res)
}
