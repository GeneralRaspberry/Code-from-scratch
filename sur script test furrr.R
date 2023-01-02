
library("purrr")
library("dplyr")



add_iteration_counter <- function(x, i) {
  # Add the iteration counter as a new column to the element x
  x$iteration <- i
  # Return the modified element
  return(x)
}

final_output_with_iteration <- map2(final_output[[2]], seq_along(final_output[[2]]), 
                                    add_iteration_counter)
#add for dataframe
final_output_df <- bind_rows(final_output_with_iteration)
timemax<-max(final_output_df$time)

#################################################################
# Define a function to combine the infected and uninfected hosts for a single landscape
combine_landscape <- function(landscape, infected) {
  # Convert the landscape to a data frame
  landscape_df <- as.data.frame(landscape)
  # Identify the uninfected hosts
  uninfected <- anti_join(landscape_df, infected, by = c("x", "y"))
  # Combine the infected and uninfected hosts
  combined <- bind_rows(infected, uninfected)
  return(combined)
}

# Apply the combine_landscape() function to each element of the landscapes list
combined_landscapes <- map(landscapes, 
                           combine_landscape, 
                           infected = final_output_df)
#This will create a list of data frames, where each data frame contains both the infected and uninfected hosts for a single landscape. You can then use bind_rows() to concatenate all of the data frames into a single data frame if desired.

combined_landscapes_df<-bind_rows(combined_landscapes)
########################attempt v2##########################################
combine_landscape <- function(landscape, infected, iteration) {
  # Convert the landscape to a data frame
  landscape_df <- as.data.frame(landscape)
  
  # Filter the infected data frame to only include the infected hosts for the current iteration
  infected_filtered <- infected[infected$iteration == iteration, ]
  
  # Identify the uninfected hosts
  uninfected <- anti_join(landscape_df, infected_filtered, by = c("x", "y"))
  
  # Assign the current iteration value to the uninfected hosts
  uninfected$iteration <- iteration
  
  # Combine the infected and uninfected hosts
  combined <- bind_rows(infected_filtered, uninfected)
  return(combined)
}

combined_landscapes <- map2(landscapes, 
                            seq_along(landscapes),
                            combine_landscape, 
                            infected = final_output_df)

combined_landscapes_df <- bind_rows(combined_landscapes)

combined_landscapes_df %>%
  filter(iteration ==4) %>%
  ggplot(aes(x = x, y = y)) +
  geom_point()



############################################################################

#surveillance script
n<-seq(15,105,length.out = 2)
fre<-seq(15,105,length.out = 2)
datalist<-list()
indicator<-1
for (f in n){
  for (l in fre){
test4<-function(x,y){
  for(g in y){
    r5<-pluck(x)
    d<-r5[sample(nrow(r5),size=f,replace=FALSE),]
    print(d)
    print(match(g,y))
print(paste0("fre=",f,"n=",l))
    
    if(g>min(d$time)){
      m<-sum(d <= g)
      q<-mean(r5$time <= g)
      #t<-min(d$time)
      #d<-match(g,y)
      mylist<-list(q,m,timedet=g,theta=final_output[[1]][[2]],beta=final_output[[1]][[1]],sim=iteration,rlf=final_output[[1]][[3]],frequency=l,samplesize=f)
      return(mylist)
    }
  }
}




stti<-sample(1:l,length(final_output_with_iteration),replace=TRUE)
samp.time <- lapply(stti, function(x) seq(from = x, to = timemax, by = l))


dftest4<-map2(final_output_with_iteration,samp.time,test4)
print("....................simulation set completed.............................")
dftestlistdocall<-data.frame(do.call(rbind,dftest4))
datalist[[indicator]]<-dftestlistdocall
      indicator<-indicator+1
  }
}

dfcompleteheatmap<-do.call(rbind,datalist)



dfcompleteheatmap$q<-as.numeric(dfcompleteheatmap$V1)
sum.q<-dfcompleteheatmap%>%group_by(beta,theta,frequency,samplesize)%>%summarise_at(vars(q),list(q_mean = mean))
sum.q<-merge(sum.q,rdataset1,by=c("theta","beta"))
sum.q$anq<-((sum.q$r_mean)*(as.numeric(sum.q$frequency)/as.numeric(sum.q$samplesize)))
sum.q$absdif<-abs(sum.q$anq-sum.q$q_mean)
sum.q$reldif<-sum.q$absdif/sum.q$q_mean
dfcompleteheatmap$t1<-as.numeric(dfcompleteheatmap$timedet)
#dfcompleteheatmap$d1<-as.numeric(dfcompleteheatmap$d)
time314<-dfcompleteheatmap%>%group_by(beta,theta,frequency,samplesize)%>%summarise_at(vars(t1),list(t_mean = mean))
#steps314<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(d1),list(steps_mean = mean))


save(dfcompleteheatmap,file="dfheatmapregsurtest1.rda")
save(sum.q,file="sum.qcompleteregsurtest1.rda")
save(rdataset1,file="rdataset1regsurtest1.Rda")
save(time314,file="time314regsurtest1.Rda")


absdifgg<-ggplot(sum.q,aes(x=as.numeric(theta),y=as.numeric(beta),fill=absdif)) +
  geom_tile() + 
  scale_fill_gradient(low="lightblue",
                      high="darkblue",
                      name="Absolute Difference")+
  theme_minimal()+
  ylab("\U03B2")+
  xlab("\U03B8")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(absdif, 3)),colour="white")

#ggsave(file="absdifregbeta7theta7.pdf")



anqgg<-ggplot(sum.q,aes(x=as.numeric(theta),y=as.numeric(beta),fill=anq)) +
  geom_tile() + 
  scale_fill_gradient(low="#FF7F7F",
                      high="darkred",
                      name="rule of thumb prediction")+
  theme_minimal()+
  ylab("\U03B2")+
  xlab("\U03B8")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(anq, 3)),colour="white")

#ggsave(file="anqvaluesregbeta7theta7.pdf")


q_meangg<-ggplot(sum.q,aes(x=as.numeric(theta),y=as.numeric(beta),fill=q_mean)) +
  geom_tile() + 
  scale_fill_gradient(low="lightgreen",
                      high="darkgreen",
                      name="Simulated detection")+
  theme_minimal()+
  ylab("\U03B2")+
  xlab("\U03B8")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(q_mean, 3)),colour="blue")

#ggsave(file="simulatedvaluesregregbeta7theta7.pdf")



reldifgg<-ggplot(sum.q,aes(x=as.numeric(theta),y=as.numeric(beta),fill=reldif)) +
  geom_tile() + 
  scale_fill_gradient(low="lightyellow",
                      high="#CCCC00",
                      name="Relative Difference")+
  theme_minimal()+
  ylab("\U03B2")+
  xlab("\U03B8")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1) + 
  geom_text(aes(label = round(reldif, 3)),colour="black")

#ggsave(file="reldifclusterregregbeta7theta7.pdf")


collectedgraphs<-ggarrange(anqgg,q_meangg,absdifgg,reldifgg,labels = c("a","b","c","d"),ncol=2,nrow=2)
ggsave(collectedgraphs,width=35,height=35,units = "cm",file="testcollectbeta50theta50surset.png")


############################time test#################################################################
dfcompleteheatmap$timedet1<-as.numeric(dfcompleteheatmap$timedet)

time364<-data.frame(dfcompleteheatmap%>%group_by(beta,theta,frequency,samplesize)%>%summarise_at(vars(timedet1),list(t_mean = mean)))
sum.q$timepred<-log(sum.q$anq/(sum.q$int))/rdataset1$r_mean
sum.q$absdif1<-abs(sum.q$timepred-time364$t_mean)
sum.q$reldif1<-sum.q$absdif1/time364$t_mean
######################################################################################################


ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=reldif1)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative time",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1) + 
  geom_text(aes(label = round(reldif1, 3)))

ggsave(file="reldif1regsurbeta150theta20.pdf")




ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=absdif1)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="abs time",
                      limits=c(0,50))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(absdif1, 3)))

ggsave(file="absdif1regsurbeta150theta20.pdf")


################################################################################
sum.q$int<-sum.q$anq/exp(rdataset1$r_mean*time364$t_mean)
sum.q$absdif2<-abs((infbegin/900)-sum.q$int)
sum.q$reldif2<-sum.q$absdif2/(infbegin/900)




ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=reldif2)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(reldif2, 3)))

ggsave(file="heatmaprel2beta150theta20start10.pdf")

ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=absdif2)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Absolute Difference initial inoculum",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1) + 
  geom_text(aes(label = round(absdif2, 3)))

ggsave(file="heatmapabs2beta150theta20start10.pdf")

###################################################################################################################

sum.q$anq1<-(log(sum.q$anq/sum.q$int)*as.numeric(sum.q$frequency))/(time364$t_mean*as.numeric(sum.q$samplesize))
sum.q$absdif3<-abs(sum.q$q_mean-sum.q$anq1)
sum.q$reldif3<-sum.q$absdif3/sum.q$q_mean

ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=absdif3)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="lunchtest",
                      limits=c(0,.25))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(absdif3, 3)))
ggsave(file="heatmapabs3beta150theta20start1.pdf")
           
ggsave(file="heatmaprel3beta150theta20start1.pdf")
##################################################################################################################

