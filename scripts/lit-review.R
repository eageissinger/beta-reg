# ---- Beta Regression Lit Review ------

# load pacakges
library(tidyverse)


lit.full<-read.delim("lit_review.txt")

head(lit.full)

lit<-lit.full%>%
  filter(!is.na(analysis.used))%>%
  select(-X)
head(lit)

# save revised file
write.csv(lit,"./lit-review-inprogress.csv",row.names = FALSE)


# load prepped data
lit<-read.csv("./lit-review/lit-review-complete.csv")%>%
  filter(method!="")# remove blank entries

lit%>%distinct(method)

# plot data
ggplot(lit)+
  geom_bar(aes(y=method),stat = "count")

# by year
ggplot(lit)+
  geom_vline(aes(xintercept=2012),linetype='dashed')+
  geom_bar(aes(x=year.published,fill=method),position = "stack",stat="count")+
  scale_fill_brewer(palette = "Set3")+
  geom_text(aes(x=year.published, fill=method, label = ..count..),stat="count",position = "stack",vjust=1.5)
# set up bar plot to be side by side, not stacked
ggplot(lit)+
  geom_bar(aes(y=year.published,fill=method),position = "fill",stat="count")+
  scale_fill_brewer(palette = "Set3")+
  geom_text(aes(y=year.published, fill=method,label = ..count..), stat="count",position="fill",hjust=2)

# summary stats
lit%>%
  group_by(method)%>%
  summarise(N=n(),timeline.min=min(year.published),timeline.max=max(year.published))%>%
  ungroup()%>%
  arrange(desc(N))

#  2012 and later (2 years after betareg package developed)
lit%>%
  filter(year.published>=2012)%>%
  group_by(method)%>%
  summarise(N=n(),timeline.min=min(year.published),timeline.max=max(year.published))%>%
  ungroup()%>%
  arrange(desc(N))

#before 2012
lit%>%
  filter(year.published<2012)%>%
  group_by(method)%>%
  summarise(N=n(),timeline.min=min(year.published),timeline.max=max(year.published))%>%
  ungroup()%>%
  arrange(desc(N))


# by discipline
lit%>%distinct(Field.Discipline)
lit%>%distinct(Tidy.Discipline)

ggplot(lit)+
  geom_bar(aes(x=year.published,fill=Tidy.Discipline),postition="stack", stat="count")+
  geom_text(aes(x=year.published, fill=Tidy.Discipline, label = ..count..),stat="count",position = "stack",vjust=1.5)+
  scale_fill_brewer(palette = "Set3")

ggplot(lit)+
  geom_vline(aes(xintercept=2012),linetype='dashed')+
  geom_bar(aes(x=year.published,fill=method),position = "stack",stat="count")+
  scale_fill_brewer(palette = "Set3")+
  #geom_text(aes(x=year.published, fill=method, label = ..count..),stat="count",position = "stack",vjust=1.5)+
  facet_wrap(~Tidy.Discipline)

lit%>%
  group_by(Tidy.Discipline)%>%
  summarise(n())
# ---- Manuscript Figures ------

#Table S1
lit%>%
  select(?..Citation,year.published,analysis.used)

# Fig 1
#png("figures/beta_carboynl-carbon_diagnostics.png",  width = 160, height = 160, units = "mm",res = 600)
Fig1<-ggplot(lit)+
  geom_bar(aes(x=year.published,fill=Tidy.Discipline),position="stack", stat="count")+
  geom_text(aes(x=year.published, fill=Tidy.Discipline, label = ..count..),stat="count",position = "stack",vjust=1.5)+
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand=c(0,0))+
  labs(y="Number of Publications",
       x="Publication Year",
       fill="Discipline")+
  theme(axis.text = element_text(size=12),
        axis.title.x=element_text(size=14,margin=margin(t=10)),
        axis.title.y=element_text(size=14,margin=margin(r=10)))+
  theme(legend.title = element_text(size=10,hjust=.5),
        legend.text = element_text(size=10),
        legend.key.size = unit(1,"cm"),
        legend.position = "top")

ggsave(plot=Fig1,filename="C:/Users/emili/Documents/Research/Beta/figures/lit-review-Fig1.png",width=160,height=160,units="mm", dpi=600)

# Fig 2
Fig2<-ggplot(lit)+
  geom_vline(aes(xintercept=2012),linetype='dashed',size=1)+
  geom_bar(aes(x=year.published,fill=method),position = "stack",stat="count")+
  scale_fill_brewer(palette = "Set3")+
  geom_text(aes(x=year.published, fill=method, label = ..count..),stat="count",position = "stack",vjust=1.5)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand=c(0,0))+
  labs(y="Number of Publications",
       x="Publication Year",
       fill="Discipline")+
  theme(axis.text = element_text(size=12),
        axis.title.x=element_text(size=14,margin=margin(t=10)),
        axis.title.y=element_text(size=14,margin=margin(r=10)))+
  theme(legend.title = element_text(size=12,hjust=.5),
        legend.text = element_text(size=10),
        legend.key.size = unit(.75,"cm"))

ggsave(plot=Fig2,filename="C:/Users/emili/Documents/Research/Beta/figures/lit-review-Fig2.png",width=180,height=160,units="mm", dpi=600)

# Table 1a
#before 2012
lit%>%
  filter(year.published<2012)%>%
  group_by(method)%>%
  summarise(N=n(),timeline.min=min(year.published),timeline.max=max(year.published))%>%
  ungroup()%>%
  arrange(desc(N))

#Table 1b
lit%>%
  filter(year.published>=2012)%>%
  group_by(method)%>%
  summarise(N=n(),timeline.min=min(year.published),timeline.max=max(year.published))%>%
  ungroup()%>%
  arrange(desc(N))
