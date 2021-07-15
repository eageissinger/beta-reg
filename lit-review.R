# ---- Beta Regression Lit Review ------

# load pacakges
library(tidyverse)

# set working directory
setwd("C:/Users/emili/Documents/Research/Beta/lit-review/")

lit.full<-read.delim("lit_review.txt")

head(lit.full)

lit<-lit.full%>%
  filter(!is.na(analysis.used))%>%
  select(-X)
head(lit)

# save revised file
write.csv(lit,"./lit-review-inprogress.csv",row.names = FALSE)


# load prepped data
lit<-read.csv("./lit-review-complete.csv")%>%
  filter(method!="")# remove blank entries

lit%>%distinct(method)

# plot data
ggplot(lit)+
  geom_bar(aes(y=method),stat = "count")

# by year
ggplot(lit)+
  geom_bar(aes(x=year.published,fill=method))
# set up bar plot to be side by side, not stacked

# summary stats
lit%>%
  group_by(method)%>%
  summarise(n(),timeline.min=min(year.published),timeline.max=max(year.published))

# by discipline
