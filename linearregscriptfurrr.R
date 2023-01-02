
library("furrr")

###############################################

newvar<-function(x){
  dftemp<-pluck(x)
  tempfunccall <- filter(dftemp, infected <= (hosts/4))#<------sync here
}

list_hosts_quartile<-map(final_output_infcount,newvar)

###############################################

# Set up the cluster
cores <- availableCores()-2
plan(multisession, workers=cores)

lineareglist <- list()
g <- 1
linearpluck <- function(x) {
  dftemp <- pluck(x)
  for (i in unique(dftemp$iteration)) {
    lmoutput <- lm(formula = log(infected) ~ time, data = filter(dftemp, iteration == i))
    lmoutput_int <- as.numeric(exp(lmoutput$coefficients[1]))
    lmoutput_r <- as.numeric(lmoutput$coefficients[2])
    lmoutput_cof <- as.numeric(summary(lmoutput)$r.squared)
    dfnew <- data.frame(int = lmoutput_int, rsquared = lmoutput_cof, r = lmoutput_r, 
                        beta = be, theta = th,
                        landscape = rlf, iteration = i)
    lineareglist[[g]] <- dfnew
    g <- g + 1
  }
  return(lineareglist)
}


  rdata <- future_map(list_hosts_quartile, linearpluck)


rdataset <- bind_rows(rdata)

#save(rdataset,file="rdatasetbeta7theta7.Rda")

##################################################################################################

newvar2<-function(mydata){
  rdatasettemp<-(pluck(mydata))
  rdatasettemp1<- unlist(lapply(rdatasettemp,`[`, c("r")))
  rdatamean<-mean(rdatasettemp1)
  rdatasd<-sd(rdatasettemp1)
  rdataframesdmean<-data.frame(r_mean=rdatamean,rdatasd=rdatasd,
                               beta=rdatasettemp[[1]]["beta"],
                               theta=rdatasettemp[[1]]["theta"],
                               landscape=rdatasettemp[[1]]["landscape"])
}


r_lnreg1<-function(i=NULL){

rdatamean<-map(rdataset,newvar2)
rdatameantable<-do.call(rbind.data.frame,rdatamean)

}
cl <- makeCluster(mc <- getOption("cl.cores", 20))
clusterCall(cl,function() library("purrr"))
clusterCall(cl,function() library("dplyr"))
clusterExport(cl=cl, varlist=c("rdataset","lineareglist","linearpluck","newvar2"),envir = environment())
par_r1<-parLapply(1,fun=r_lnreg1,cl=cl)
stopCluster(cl)
rdataset1 <- do.call("rbind", par_r1)




#rdatamean<-rdataset%>%group_by(beta, theta,landscape)%>%summarise_at(vars(r),list(r_mean = mean))
#library(tidyverse)
#tempLongPreds <- pivot_longer(temp, cols=c(infected,predExp,predLog), names_to="estimate", values_to="inf")
ggplot(rdataset1)+geom_line(aes(x=theta,y=r_mean,group=as.factor(beta),colour=as.factor(beta)))
ggplotrexp1<-ggplot(rdataset1)+geom_line(aes(x=theta,y=r_mean,group=as.factor(beta),colour=as.factor(beta)))
