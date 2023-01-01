library("purrr")
library("dplyr")

#function to count how many hosts were infected at each time stamp
infcount<-function(x,beta, theta, landscape,i){
  datasp<-pluck(x)
  times <- sort(unique(datasp$time))
  infcountframes<-datasp  %>% 
    do(data.frame(iteration=i,beta=beta,theta=theta,landscape=landscape,
                  time=times, infected=sapply(times, function(x) sum(.$time <= x))))
}
#we use map2 to mutate() as anonymous the iteration column
final_output_infcount <- map2(final_output[[2]], seq_along(final_output[[2]]), 
                              infcount, 
                              beta = final_output[[1]][[1]], 
                              theta = final_output[[1]][[2]], 
                              landscape = final_output[[1]][[3]])
#dataframe of results
final_output_infcount_df <- bind_rows(final_output_infcount)


# Add a counter column to each element in the data frame to track simulation
final_output_infcount_df <- final_output_infcount_df %>%
  group_by(list_element) %>%
  mutate(counter = row_number())