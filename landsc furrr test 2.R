rm(list=ls)
library("spatstat")
library("furrr")
library("tidyverse")
set.seed(seed=NULL)

numsims<-1000

episet<-as.list(seq(1:numsims))

generate_landscape<-function(radiusCluster,
                            lambdaParent,
                             lambdaDaughter,
                             randmod,
                             hosts,
                             dim){
numbparents<-rpois(1,lambdaParent*dim)

xxParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
yyParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)

numbdaughter<-rpois(numbparents,(lambdaDaughter))
sumdaughter<-sum(numbdaughter)



thetaLandscape<-2*pi*runif(sumdaughter)

rho<-radiusCluster*sqrt(runif(sumdaughter))



xx0=rho*cos(thetaLandscape)
yy0=rho*sin(thetaLandscape)


xx<-rep(xxParent,numbdaughter)
yy<-rep(yyParent,numbdaughter)

xx<-xx+xx0

yy<-yy+yy0
cds<-data.frame(xx,yy)
is_outlier<-function(x){
  x > dim| x < 0
}
cds<-cds[!(is_outlier(cds$xx)|is_outlier(cds$yy)),]
while (nrow(cds)<hosts){
  dif<-hosts-nrow(cds)
  extraparentxx<-sample(xxParent,dif,replace = TRUE)
  extraparentyy<-sample(yyParent,dif,replace = TRUE)
  extrathetaLandscape<-2*pi*runif(dif)
  extrarho<-radiusCluster*sqrt(runif(dif))
  newextracoodsxx<-extrarho*cos(extrathetaLandscape)
  newextracoodsyy<-extrarho*sin(extrathetaLandscape)
  extraxx<-extraparentxx+newextracoodsxx
  extrayy<-extraparentyy+newextracoodsyy
  cdsextra<-data.frame(xx=extraxx,yy=extrayy)
  cds<-rbind(cds,cdsextra)
}

sampleselect<-sample(1:nrow(cds),hosts,replace=F)
cds<-cds%>%slice(sampleselect)

randfunction<-function(x){
  x<-runif(length(x),0,dim)
}
randselect<-sample(1:nrow(cds),floor(hosts*randmod),replace=F)
cds[randselect,]<-apply(cds[randselect,],1,randfunction)

landscape<-ppp(x=cds$xx,y=cds$yy,window=owin(xrange=c(0,dim),yrange=c(0,dim)))
landscapestore<-data.frame(landscape)
landscapestore
}

# Use 4 cores to generate the landscapes in parallel
cores <- availableCores()-2
plan(multisession, workers=cores)

# Generate the landscapes in parallel
landscapes <- future_map(episet, generate_landscape,
                         lambdaParent=.05,
                         lambdaDaughter=25,
                         randmod=0,
                         hosts=1000,
                         dim=1000)
#rlang::last_error()
#plot(landscape)
#print(length(landscape$marks))
