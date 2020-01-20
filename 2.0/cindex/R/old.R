old_cindex=function(yhat,y,status,w=rep(1,length(y))){
  # Rob's implementation of C-index
  #note:  this function gives identical results to the summary function from coxph, and the concordance.index function in survcomp.
  #  (with their default settings), and no ties in yhat. Works with ties in y. But does not agree with latter when yhat has ties. There are conflicting definitions for c-index in this case
  #
  #formula used  is
  #  Concordance = (#all concordant pairs + #tied pairs/2)/(#total pairs including ties).
  # with w,  weights used are ave wts for each pair
  time.computeC.start <- Sys.time()

  risksets=which(status==1)
  w=length(w)*w/sum(w)
  fun=function(riskset,y,yhat,w){
    total=concordant=0
    i=riskset
    rest=which(y>y[i])
    if(length(rest)>0){
      ww=w[rest]*w[i]
      total=sum(ww)
      concordant = 1.0*sum(ww*(yhat[rest]<yhat[i]))+0.5*sum(ww*(yhat[rest]==yhat[i]))
    }
    return(c(concordant,total))
  }
  out=sapply(risksets,fun,y,yhat,w)
  cindex=sum(out[1,])/sum(out[2,])

  return(cindex)
}
