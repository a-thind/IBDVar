library(dplyr)
library(ggplot2)
library(RColorBrewer)
# summary functions
counts <- function(data, feature) {
  data %>%
    group_by({{ feature }}) %>%
    summarise(count = n()) %>% arrange(desc(count))
}

# produce barplot showing counts for each type of group variable
plot_counts <- function(data, var_name, flip=FALSE, title=NULL) {
  if (flip) {
    ggplot(data, aes(x=pull(data, 1), y=count, pull(data, 1))) +
      geom_bar(stat="identity") +
      xlab(var_name) +
      coord_flip() +
      ylab("Count") +
      scale_fill_discrete(name=var_name) + 
      ggtitle(title)
  } else {
    ggplot(data, aes(x=pull(data, 1), y=count, fill=pull(data, 1))) +
      geom_bar(stat="identity") +
      xlab(var_name) +
      ylab("Count") +
      ggtitle(title) +
      scale_fill_discrete(name=var_name)
  }
  
}