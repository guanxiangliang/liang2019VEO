#### Figure 2A ####
mtmp<-dna.rpm.family%>%dplyr::select(library_id,uniprot.family,rpm)
names(mtmp)<-c("library_id","species","rpm")
df<-merge(mtmp,meta)
tmp<-df%>%filter(study_group!="Negative control",species=="Ackermannviridae"|species=="Anelloviridae"|species=="Herelleviridae"|species=="Inoviridae"|species=="Mimiviridae"|species=="Podoviridae"|species=="Myoviridae"|species=="Siphoviridae"|species=="Microviridae")
tmp$species<-factor(tmp$species,levels=c("Ackermannviridae","Anelloviridae","Herelleviridae","Inoviridae","Mimiviridae","Podoviridae","Myoviridae","Siphoviridae","Microviridae"))
ggplot(tmp,aes(species,rpm/1000000,fill=species))+
  geom_boxplot(outlier.size = 0.5)+
  #geom_beeswarm(cex=1,alpha=0.8,size=1,shape=21,colour="white",fill="black",groupOnX = T)+
  labs(x=c(""),y=c("Abundance"))+
  theme_classic()+
  coord_flip()+
  theme(axis.text.y = element_text(size=20,colour = "black"),
        axis.title = element_text(size=25,colour="black"),
        axis.text.x = element_text(size=20,colour = "black"),
        legend.position = "none")
ggsave("output/main.figure/DNA.taxa.pdf",width=7,height =5)
#### End ####

#### Extra analysis for comments ####
tmp1<-df%>%dplyr::select(library_id,rpm,study_group,species)%>%
  filter(study_group!="NC",species=="Ackermannviridae"|species=="Herelleviridae"|species=="Podoviridae"|species=="Myoviridae"|species=="Siphoviridae")%>%
  group_by(library_id) %>%
  summarise(abun=sum(rpm))%>%
  ungroup() %>%
  mutate(family="Caudovirales")
tmp2<-df%>%dplyr::select(library_id,rpm,study_group,species)%>%
  filter(study_group!="NC",species=="Microviridae")%>%
  group_by(library_id) %>%
  summarise(abun=sum(rpm))%>%
  ungroup()%>%
  mutate(family="Microviridae")
tmp<-rbind (tmp1,tmp2) %>%inner_join(meta)
df<-tmp%>%mutate(presence=ifelse(abun>10000,"Yes","No"))
df<-df %>% filter(study_group=="VEO"|study_group=="Healthy",family=="Microviridae") %>% select(study_group,presence)
df$study_group<-droplevels(df$study_group)
fisher.test(table(df))[1][[1]]

tmp<-rbind (tmp1,tmp2) %>%inner_join(meta)
df<-tmp%>%mutate(presence=ifelse(abun>10000,"Yes","No"))
df<-df %>% filter(study_group=="VEO"|study_group=="Healthy",family=="Caudovirales") %>% select(study_group,presence)
df$study_group<-droplevels(df$study_group)
fisher.test(table(df))[1][[1]]
#### End ####

#### Figure 2B, 2C and 2D ####
tmp1<-df%>%dplyr::select(library_id,rpm,study_group,species)%>%
  filter(study_group!="NC",species=="Ackermannviridae"|species=="Herelleviridae"|species=="Podoviridae"|species=="Myoviridae"|species=="Siphoviridae")%>%
  group_by(library_id) %>%
  summarise(abun=sum(rpm))%>%
  ungroup() %>%
  mutate(family="Caudovirales")
tmp2<-df%>%dplyr::select(library_id,rpm,study_group,species)%>%
  filter(study_group!="NC",species=="Microviridae")%>%
  group_by(library_id) %>%
  summarise(abun=sum(rpm))%>%
  ungroup()%>%
  mutate(family="Microviridae")
tmp<-rbind (tmp1,tmp2) %>%inner_join(meta)

# 1C
tmp$study_group<-factor(tmp$study_group,levels=c("VEO-IBD","Healthy controls"))

tmp$cal_group<-factor(tmp$cal_group,levels=c("Active VEO-IBD","Inactive VEO-IBD"))

p<-ggplot(tmp,aes(x=family,y=abun/1000000,color=family))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_beeswarm(data=.%>%filter(family=="Caudovirales"),cex=2,alpha=0.8,size=3,shape=21,colour="white",fill=guan.color[1],groupOnX = T)+
  geom_beeswarm(data=.%>%filter(family=="Microviridae"),cex=2,alpha=0.8,size=3,shape=21,colour="white",fill=guan.color[2],groupOnX = T)+
  scale_color_manual(values = guan.color[1:2])+
  expand_limits(y = c(0, 1.1))+
  labs(x=c(""),y=c("Abundance"))+
  facet_grid(.~study_group)+
  theme_classic()+
  theme(axis.text = element_text(size=20,colour = "black"),
        axis.title = element_text(size=25,colour="black"),
        strip.text = element_text(size=30,colour="black"),
        legend.position = "none")

g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
fills <- guan.color[c(3,5)]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


pdf("output/main.figure/Caudo.vs.Micro.pdf", width = 10, height = 5)
grid.draw(g)
grid.text(c("p = 0.03","p = 0.34"),gp=gpar(fontsize=15), x = c(0.31,0.77), y = c(0.85,0.85))
grid.segments(x0=c(0.22,0.68),x1=c(0.40,0.86),y0=c(0.81,0.81),y1=c(0.81,0.81))
dev.off()

# 1D


p<-ggplot(tmp[tmp$study_group=="VEO-IBD"&!is.na(tmp$cal_group),],aes(x=family,y=abun/1000000,color=family))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_beeswarm(data=.%>%filter(family=="Caudovirales"),cex=2,alpha=0.8,size=3,shape=21,colour="white",fill=guan.color[1],groupOnX = T)+
  geom_beeswarm(data=.%>%filter(family=="Microviridae"),cex=2,alpha=0.8,size=3,shape=21,colour="white",fill=guan.color[2],groupOnX = T)+
  scale_color_manual(values = guan.color[1:2])+
  expand_limits(y = c(0, 1.1))+
  labs(x=c(""),y=c("Abundance"))+
  facet_grid(.~cal_group)+
  theme_classic()+
  theme(axis.text = element_text(size=20,colour = "black"),
        axis.title = element_text(size=25,colour="black"),
        strip.text = element_text(size=30,colour="black"),
        legend.position = "none")

g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
fills <- guan.color[c(4,6)]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/main.figure/Caudo.vs.Micro.cal_group.pdf", width = 10, height = 5)
grid.draw(g)
grid.text(c("p = 0.02","p = 0.96"),gp=gpar(fontsize=15), x = c(0.31,0.77), y = c(0.85,0.85))
grid.segments(x0=c(0.22,0.68),x1=c(0.40,0.86),y0=c(0.81,0.81),y1=c(0.81,0.81))
dev.off()

# gender
wilcox.test(tmp[tmp$Gender=="Female"&tmp$family=="Caudovirales",]$abun,tmp[tmp$Gender=="Female"&tmp$family=="Microviridae",]$abun)
wilcox.test(tmp[tmp$Gender=="Male"&tmp$family=="Caudovirales",]$abun,tmp[tmp$Gender=="Male"&tmp$family=="Microviridae",]$abun)

# compare study group
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae",]$abun)
wilcox.test(tmp[tmp$study_group=="Healthy controls"&tmp$family=="Caudovirales",]$abun,tmp[tmp$study_group=="Healthy controls"&tmp$family=="Microviridae",]$abun)
# disease activity
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$cal_group=="Active VEO-IBD"&tmp$family=="Caudovirales",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$cal_group=="Active VEO-IBD"&tmp$family=="Microviridae",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$cal_group=="Inactive VEO-IBD"&tmp$family=="Caudovirales",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$cal_group=="Inactive VEO-IBD"&tmp$family=="Microviridae",]$abun)
# immunosuppression
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="Treatment",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="Treatment",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="No treatment",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="No treatment",]$abun)
# un immunosuppression group compare activity
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="No treatment"&tmp$cal_group=="Active VEO-IBD",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="No treatment"&tmp$cal_group=="Active VEO-IBD",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="No treatment"&tmp$cal_group=="Inactive VEO-IBD",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="No treatment"&tmp$cal_group=="Inactive VEO-IBD",]$abun)
# inactive group compare immunosuppression
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="Yes"&tmp$cal_group=="Inactive VEO-IBD",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="Yes"&tmp$cal_group=="Inactive VEO-IBD",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="No"&tmp$cal_group=="Inactive VEO-IBD",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="No"&tmp$cal_group=="Inactive VEO-IBD",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="Yes"&tmp$cal_group=="Active VEO-IBD",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="Yes"&tmp$cal_group=="Active VEO-IBD",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$immu_suppress=="No"&tmp$cal_group=="Active VEO-IBD",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$immu_suppress=="No"&tmp$cal_group=="Active VEO-IBD",]$abun)
# antibiotics
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$current_antibiotics=="Treatment",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$current_antibiotics=="Treatment",]$abun)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales"&tmp$current_antibiotics=="No treatment",]$abun,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae"&tmp$current_antibiotics=="No treatment",]$abun)



p<-ggplot(tmp[tmp$study_group=="VEO-IBD"&!is.na(tmp$cal_group),],aes(x=family,y=abun/1000000,color=family))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_beeswarm(data=.%>%filter(family=="Caudovirales"),cex=2,alpha=0.8,size=3,shape=21,colour="white",fill=guan.color[1],groupOnX = T)+
  geom_beeswarm(data=.%>%filter(family=="Microviridae"),cex=2,alpha=0.8,size=3,shape=21,colour="white",fill=guan.color[2],groupOnX = T)+
  scale_color_manual(values = guan.color[1:2])+
  expand_limits(y = c(0, 1.1))+
  labs(x=c(""),y=c("Relative abundance"))+
  facet_grid(.~immu_suppress)+
  theme_classic()+
  theme(axis.text = element_text(size=20,colour = "black"),
        axis.title = element_text(size=25,colour="black"),
        strip.text = element_text(size=30,colour="black"),
        legend.position = "none")

g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
fills <- guan.color[c(9,10)]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}


pdf("/media/lorax/users/guanxiang/2_veo/liang.veo/output/Caudo.vs.Micro.immunosuppression.pdf", width = 10, height = 5)
grid.draw(g)
grid.text(c("p = 0.002","p = 0.94"),gp=gpar(fontsize=15), x = c(0.31,0.77), y = c(0.85,0.85))
grid.segments(x0=c(0.22,0.68),x1=c(0.40,0.86),y0=c(0.81,0.81),y1=c(0.81,0.81))
dev.off()

summary(lm(abun~age,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Caudovirales",]))
summary(lm(abun~age,tmp[tmp$study_group=="Healthy controls"&tmp$family=="Caudovirales",]))

summary(lm(abun~age,tmp[tmp$study_group=="VEO-IBD"&tmp$family=="Microviridae",]))
summary(lm(abun~age,tmp[tmp$study_group=="Healthy controls"&tmp$family=="Microviridae",]))

# 1B
p<- tmp%>%filter(study_group=="Healthy controls"|study_group=="VEO-IBD")%>%
  dplyr::select(library_id,abun,family,study_group) %>%
  spread(key=family,value=abun)%>%
  ggplot(aes(x=Caudovirales/1000000,y=Microviridae/1000000,color=study_group))+
  geom_point(size=3)+
  geom_smooth(aes(x=Caudovirales/1000000,y=Microviridae/1000000),inherit.aes = F,method = "lm",colour="black",linetype = "dashed")+
  scale_color_manual(values = guan.color[c(3,5)])+
  labs(x=c("Caudovirales"),y=c("Microviride"))+
  annotate("text", x=0.75, y=0.65,size=6, label= c("p < 2.2e-16"))+
  annotate("text", x=0.75, y=0.58,size=6, label= c("r = -0.80"),hjust = 0.65)+
  theme_classic()+
  theme(axis.text = element_text(size=20,colour = "black"),
        axis.title = element_text(size=25,colour="black"),
        legend.position = c(0.75,0.8),
        legend.title = element_blank(),
        legend.text = element_text(size=15,colour = "black")
        )
pdf("output/main.figure/Caudo.vs.Micro.correlation.pdf", width = 5, height = 5,useDingbats = F)
p
dev.off()



# ratio test
tmp<-inner_join(tmp1%>%rename(Caudovirales=abun)%>%select(-family),tmp2%>%rename(Microviridae=abun)%>%select(-family))%>%
  mutate(ratio=Caudovirales/Microviridae)%>%
  select(library_id,ratio)%>%
  inner_join(meta)
tmp$study_group<-gsub(tmp$study_group,pattern = "VEO",replacement = "VEO-IBD")
tmp$study_group<-gsub(tmp$study_group,pattern = "Healthy",replacement = "Healthy controls")
tmp$study_group<-factor(tmp$study_group,levels=c("VEO-IBD","Healthy controls"))


wilcox.test(tmp[tmp$study_group=="VEO-IBD",]$ratio,tmp[tmp$study_group=="Healthy controls",]$ratio)
wilcox.test(tmp[tmp$Gender=="Female",]$ratio,tmp[tmp$Gender=="Male",]$ratio)
wilcox.test((tmp%>%filter(study_group=="VEO-IBD",cal_group=="Active"))$ratio,(tmp%>%filter(study_group=="VEO-IBD",cal_group=="Inactive"))$ratio)
wilcox.test((tmp%>%filter(study_group=="VEO-IBD",immu_suppress=="No",cal_group=="Active"))$ratio,(tmp%>%filter(study_group=="VEO-IBD",immu_suppress=="No",cal_group=="Inactive"))$ratio)


wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$immu_suppress=="Yes",]$ratio,tmp[tmp$study_group=="VEO-IBD"&tmp$immu_suppress=="No",]$ratio)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$anti.TNF=="Yes",]$ratio,tmp[tmp$study_group=="VEO-IBD"&tmp$anti.TNF=="No",]$ratio)
wilcox.test(tmp[tmp$study_group=="VEO-IBD"&tmp$current_antibiotics=="Yes",]$ratio,tmp[tmp$study_group=="VEO-IBD"&tmp$current_antibiotics=="No",]$ratio)

summary(lm(ratio~age,tmp[tmp$study_group=="VEO-IBD",]))
summary(lm(ratio~cal_value,tmp[tmp$study_group=="VEO-IBD",]))
summary(lm(ratio~time_diagno,tmp[tmp$study_group=="VEO-IBD",]))



ggplot(tmp%>%filter(study_group=="VEO-IBD",!is.na(cal_group)),aes(x=cal_group,y=ratio))+
  geom_boxplot(outlier.shape = NA,size=1)+
  scale_y_continuous(trans='log2')+
  geom_point()

ggplot(tmp,aes(x=study_group,y=ratio))+
  geom_boxplot(outlier.shape = NA,size=1)+
  scale_y_continuous(trans='log2')+
  geom_point()

#### End ####

#### Figure 3 ####

mtmp<-dna.rpm.species%>%select(library_id,uniprot.species,rpm)
names(mtmp)<-c("library_id","species","rpm")
df<-merge(mtmp,meta)

dim(com.study)
dim(com.cal)

de.spe<-inner_join(com.study,com.cal,by="uniprot.species")
tmp1<-df%>%filter(species%in%de.spe$uniprot.species)%>%
  group_by(species,study_group)%>%
  summarise(median=median(rpm))%>%
  spread(key=study_group,value=median)%>%
  mutate(de=ifelse(`VEO-IBD`>=`Healthy controls`,"VEO-IBD","Healthy controls"),fc=`VEO-IBD`/`Healthy controls`)

tmp2<-df%>%filter(species%in%de.spe$uniprot.species)%>%
  group_by(species,cal_group)%>%
  summarise(median=median(rpm))%>%
  spread(key=cal_group,value=median)%>%
  mutate(de=ifelse(`Active VEO-IBD`>=`Inactive VEO-IBD`,"VEO-IBD","Healthy controls"),fc=`Active VEO-IBD`/`Inactive VEO-IBD`)

de.spe<-inner_join(tmp1,tmp2,by="species")%>%mutate(select=ifelse(de.x==de.y,"Yes","No"))%>%
  filter(select=="Yes")%>%
  filter(fc.x>1.5|fc.x<0.67)%>%
  filter(fc.y>1.5|fc.y<0.67)



# 1. species study_group 
tmp<-df%>%
  filter(study_group=="VEO-IBD")%>%
  group_by(species)%>%
  summarise(median=median(rpm))%>%
  arrange(median)
df$species<-factor(df$species,levels=tmp$species,ordered = T)

tmp<-df%>%filter(species %in% as.character(de.spe$species),!is.na(cal_group))%>%
  mutate(de=ifelse(species=="Brochothrix phage BL3","VEO-IBD","Healthy controls"))

g<-ggplot(tmp,aes(x=species,y=log(rpm+1,10),fill=cal_group))+
  geom_boxplot(outlier.shape = NA)+
  #scale_y_continuous(position = "right")+
  coord_cartesian(xlim=c(0, 1e4))+
  scale_fill_manual(values=guan.color[c(4,6,5)])+
  coord_flip(ylim=c(0, 4))+
  ylab(bquote(''*Log[10]*'RPM'))+
  xlab("")+
  facet_grid(de~.,scales = "free",space = "free_y")+
  theme_classic()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title = element_text(size=15,color="black"),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=15,colour = "black")) # subscript in ggplot

pdf("output/main.figure/de.speceis.pdf", width = 10, height = 5)
g
dev.off()

#### End ####

#### Figure 4 ####
df<-dna.filter%>%filter(family=="Anelloviridae")
df<-df%>%inner_join(meta)
df$infection<-"Uninfected"
df[df$cov>=0.33,"infection"]<-"Infected"
tmp<-df %>% filter(study_group=="VEO"|study_group=="Healthy") %>% dplyr::select(study_group,infection)
tmp$study_group<-droplevels(tmp$study_group)
fisher.test(table(tmp))[1][[1]]

g<-tmp%>%
  group_by(study_group)%>%
  count(infection)%>%
  mutate(ratio=n/sum(n))%>%
  filter(infection=="Infected")%>%
ggplot() +
  #geom_errorbar(aes(ymin=ratio*100, ymax=se*100,x = Infant_feeding_type), width=0.2,size=0.5,
  #              position=position_dodge(.9)) +
  geom_bar(aes(y = ratio*100, x = study_group,fill=study_group), color="black",stat="identity",width = 0.5,size=0.5)+
  scale_fill_manual(values= guan.color[c(3,5)])+
  expand_limits(y = c(0, 32))+
  labs(title="",x ="", y = "Percentage of subjects (%)\nwith Anelloviridae")+
  scale_x_discrete(labels=c( "VEO-IBD","Healthy controls"))+
  #facet_grid(.~Cohort,scales = "free", space = "free")+
  annotate("text", x = c(1.5), y = c(30), size=5, label= c("p = 0.03"))+
  annotate("segment", x = c(1.05), xend = c(1.95), y = c(29), yend = c(29))+
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text.x = element_text(size=20,colour = "black"),
        axis.text.y = element_text(size=20,colour = "black",margin = margin(t = -20, unit = "pt")),
        strip.text = element_text(size=25),
        #strip.background = element_rect(fill=cohort.col[2]),
        legend.position = "none"
  )
pdf("output/main.figure/anelloviridae.study_group.pdf", width = 5, height = 5)
g
dev.off()




tmp<-df %>% filter(cal_group=="Yes"|cal_group=="No") %>% select(cal_group,infection)
tmp$cal_group<-droplevels(tmp$cal_group)
fisher.test(table(tmp))[1][[1]]


tmp<-df %>% filter(immu_suppress=="Yes"|immu_suppress=="No",study_group=="VEO") %>% select(immu_suppress,infection)
tmp$immu_suppress<-gsub(tmp$immu_suppress,pattern = "Yes",replacement = "Treatment")
tmp$immu_suppress<-gsub(tmp$immu_suppress,pattern = "No",replacement = "No treatment")
tmp$immu_suppress<-factor(tmp$immu_suppress,levels=c("Treatment","No treatment"))
tmp$immu_suppress<-droplevels(tmp$immu_suppress)
fisher.test(table(tmp))[1][[1]]

g<-tmp%>%
  group_by(immu_suppress)%>%
  count(infection)%>%
  mutate(ratio=n/sum(n))%>%
  filter(infection=="Infected")%>%
  ggplot() +
  #geom_errorbar(aes(ymin=ratio*100, ymax=se*100,x = Infant_feeding_type), width=0.2,size=0.5,
  #              position=position_dodge(.9)) +
  geom_bar(aes(y = ratio*100, x = immu_suppress,fill=immu_suppress), color="black",stat="identity",width = 0.5,size=0.5)+
  scale_fill_manual(values= guan.color[c(9,10)])+
  expand_limits(y = c(0,50))+
  labs(title="",x ="", y = "Percentage of subjects (%)\nwith Anelloviridae")+
  #scale_x_discrete(labels=c( "VEO-IBD","Healthy controls"))+
  #facet_grid(.~Cohort,scales = "free", space = "free")+
  annotate("text", x = c(1.5), y = c(46), size=5, label= c("p = 0.02"))+
  annotate("segment", x = c(1.05), xend = c(1.95), y = c(44.5), yend = c(44.5))+
  theme_classic()+
  theme(axis.title = element_text(size=25,colour="black"),
        axis.text.x = element_text(size=20,colour = "black"),
        axis.text.y = element_text(size=20,colour = "black",margin = margin(t = -20, unit = "pt")),
        strip.text = element_text(size=25),
        #strip.background = element_rect(fill=cohort.col[2]),
        legend.position = "none"
  )
pdf("output/main.figure/anelloviridae.treatment_group.pdf", width = 5, height = 5)
g
dev.off()

#### End ####

#### Figure 5A ####
#read<-read_delim("./input/kmer/veo.ibd/rna/veo.ibd.rna.total",delim="\t")

kmer<-read.csv("./input/kmer/veo.ibd/rna/veo.ibd.rna.pvalue.csv")

tmp<-kmer%>%filter(enrichment=="VEO"|enrichment=="Healthy")
tmp$p<-as.character(tmp$p)
tmp$p<-as.numeric(tmp$p)
tmp<-mutate(tmp,adjusted.p=p.adjust(tmp$p, method = "BH", n = length(tmp$p))) %>%
  arrange(p) %>% filter(enrichment=="VEO")
res<-tmp%>%filter(adjusted.p<0.05)  
key<-left_join(res,read,by = c("kmer" = "id"))


tmp<-meta%>%filter(study_group=="Healthy"|study_group=="VEO",library=="RNA")
tmp<-key%>%dplyr::select(as.character(tmp$library_id),kmer)%>%
  gather(key=library_id,value=reads,-kmer)%>%
  mutate(pre=ifelse(reads>0,"1","0"))%>%
  inner_join(meta)

tmp<-tmp[!duplicated(tmp[c("library_id","kmer")]),]


nt<-read.delim("input/kmer/veo.ibd/rna/veo.ibd.rna.enrichment.nt",header=F)
names(nt)<-c("kmer", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
nt<-nt[!duplicated(nt$kmer),][,c("kmer","sseqid")]
names(nt)<-c("kmer","nt.top.hit")
nt$nt.top.hit<-separate(data = nt, col = nt.top.hit, into = c("1", "2","3","ntname"),
                        sep = "\\|")$ntname
taxaId<-accessionToTaxa(nt$nt.top.hit,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql")
taxa<-getTaxonomy(taxaId,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql")
nt<-cbind(nt,taxa)

tmp<-tmp%>%left_join(nt,by="kmer")


tmp1<-tmp%>%dplyr::select(kmer,pre,library_id)%>%
  spread(key=library_id,value=pre)


tmp2<-as.data.frame(rowSums(tmp1[,-1]>0))
names(tmp2)<-"n"
tmp2$name<-tmp1[,1]

tmp$kmer <- factor(x = tmp$kmer,
                     levels = (tmp2%>%arrange(n))$name,
                     ordered = TRUE)

tmp2<-as.data.frame(colSums(tmp1[,-1]>0))
names(tmp2)<-"n"
tmp2$name<-rownames(tmp2)

tmp$library_id <- factor(x = tmp$library_id,
                         levels = (tmp2%>%arrange(desc(n)))$name,
                         ordered = TRUE)

tmp$study_group<-gsub(tmp$study_group,pattern = "VEO",replacement = "VEO-IBD")
tmp$study_group<-gsub(tmp$study_group,pattern = "Healthy",replacement = "Healthy controls")
tmp$study_group<-factor(tmp$study_group,levels=c("VEO-IBD","Healthy controls"))
tmp$pre<-factor(tmp$pre,levels=c(1,0))

tmp$superkingdom<-as.character(tmp$superkingdom)
tmp$superkingdom<-ifelse(is.na(tmp$superkingdom), "Unknown", tmp$superkingdom)
tmp$superkingdom<-factor(tmp$superkingdom,levels=c("Bacteria","Eukaryota","Unknown"))

p<-ggplot(tmp,aes(x=library_id,y=kmer,fill=pre))+
  geom_tile()+
  facet_grid(superkingdom~study_group,scales = "free",space = "free",switch = "y")+
  scale_fill_manual(values = c(guan.color[11],"white"))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "grey80",fill = NA,size = 0.1),
        #panel.margin.x=unit(0.5, "lines") , 
        #panel.border = element_blank(),
        panel.spacing.y = unit(0,"lines"),
        legend.position = "none",
        
        strip.text.x  = element_text(colour = "black",size=25),
        strip.text.y  = element_text(colour = "black",size=0))

g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
fills <- guan.color[c(3,5)]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

stripl <- which(grepl('strip-l', g$layout$name))
fills <- palette[1:5]
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/main.figure/VEO-IBD.DE.kmer.pdf", width = 10, height = 5)
grid.draw(g)
dev.off()
#### End ####

#### Figure 5B ####
#read<-read_delim("./input/kmer/veo.ibd/rna/veo.ibd.rna.total",delim="\t")

ihmp.meta<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/ihmp/ibd.virome.meta.txt")
key<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/ihmp/ibd.ihmp.enrichment.total")



tmp<-key%>%
  gather(key=library_id,value=reads,-kmer)%>%
  mutate(pre=ifelse(reads>0,"1","0"))%>%
  inner_join(ihmp.meta)


nt<-read.delim("./input/kmer/ihmp/ibd.ihmp.enrichment.nt",header=F)
names(nt)<-c("kmer", "sseqid", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
nt<-nt[!duplicated(nt$kmer),][,c("kmer","sseqid")]
names(nt)<-c("kmer","nt.top.hit")
nt$nt.top.hit<-separate(data = nt, col = nt.top.hit, into = c("1", "2","3","ntname"),
                        sep = "\\|")$ntname
taxaId<-accessionToTaxa(nt$nt.top.hit,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql")
taxa<-getTaxonomy(taxaId,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql")
nt<-cbind(nt,taxa)

tmp<-tmp%>%left_join(nt,by="kmer")
#tmp<-tmp[!duplicated(tmp[c("library_id","kmer")]),]
tmp1<-tmp%>%dplyr::select(kmer,pre,library_id)%>%
  spread(key=library_id,value=pre)


tmp2<-as.data.frame(rowSums(tmp1[,-1]>0))
names(tmp2)<-"n"
tmp2$name<-tmp1[,1]

tmp$kmer <- factor(x = tmp$kmer,
                   levels = (tmp2%>%arrange(n))$name,
                   ordered = TRUE)

tmp2<-as.data.frame(colSums(tmp1[,-1]>0))
names(tmp2)<-"n"
tmp2$name<-rownames(tmp2)

tmp$library_id <- factor(x = tmp$library_id,
                         levels = (tmp2%>%arrange(desc(n)))$name,
                         ordered = TRUE)

tmp$study_group<-gsub(tmp$study_group,pattern = "VEO",replacement = "IBD")
tmp$study_group<-gsub(tmp$study_group,pattern = "Healthy",replacement = "Healthy controls")
tmp$study_group<-factor(tmp$study_group,levels=c("IBD","Healthy controls"))
tmp$pre<-factor(tmp$pre,levels=c(1,0))

tmp$superkingdom<-as.character(tmp$superkingdom)
tmp$superkingdom<-ifelse(is.na(tmp$superkingdom), "Unknown", tmp$superkingdom)
tmp$superkingdom<-factor(tmp$superkingdom,levels=c("Bacteria","Eukaryota","Unknown"))


p<-ggplot(tmp,aes(x=library_id,y=kmer,fill=pre))+
  geom_tile()+
  facet_grid(superkingdom~study_group,scales = "free",space = "free",switch = "y")+
  scale_fill_manual(values = c(guan.color[8],"white"))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "grey80",fill = NA,size = 0.1),
        #panel.margin.x=unit(0.5, "lines") , 
        #panel.border = element_blank(),
        panel.spacing.y = unit(0,"lines"),
        legend.position = "none",
        
        strip.text.x  = element_text(colour = "black",size=25),
        strip.text.y  = element_text(colour = "black",size=0))

g <- ggplot_gtable(ggplot_build(p))
stript <- which(grepl('strip-t', g$layout$name))
fills <- guan.color[c(3,5)]
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

stripl <- which(grepl('strip-l', g$layout$name))
fills <- palette[1:5]
k <- 1
for (i in stripl) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("output/main.figure/iHMP-IBD.DE.kmer.pdf", width = 10, height = 5)
#tiff("output/main.figure/iHMP-IBD.DE.kmer.tif", compression = "lzw",width = 10, height = 5)
grid.draw(g)
dev.off()
#### End ####

#### de.id.ihmp ####
ihmp.meta<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/ihmp/ibd.virome.meta.txt")
key<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/ihmp/ibd.ihmp.enrichment.total")
de.id<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/ihmp/ibd.ihmp.de.id",header=F)

tmp<-key%>%
  filter(kmer %in% gsub(as.character(de.id$V1),replacement = "veo.ibd.rna",pattern = "ibd.ihmp")) %>%
  gather(key = library_id,value=reads,-kmer)%>%
  left_join(ihmp.meta)%>%
  group_by(study_group)%>%
  droplevels()%>%
  ggplot(aes(x=study_group,y=log(reads+1,10)))+
  geom_boxplot()+
  facet_grid(.~kmer,scales = "free_y")
  
  
  
  
  


#### End ####

#### de.id.veo.ibd ####
de.id<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/veo.ibd/rna/de.id",header=F)
#key<-read.delim("/media/lorax/users/guanxiang/2_veo/liang.veo/input/kmer/veo.ibd/rna/veo.ibd.rna.total")


#veo.ibd.key<-key%>%
#  filter(id %in% as.character(de.id$V1)) %>%
#  gather(key = library_id,value=reads,-id)%>%
#  left_join(meta)%>%
#  group_by(study_group)%>%
#  droplevels()

  ggplot(veo.ibd.key,aes(x=study_group,y=reads))+
  geom_boxplot()+
  facet_grid(.~id,scales = "free_y")
  
  ggplot(veo.ibd.key,aes(x=study_group,y=reads))+
    geom_boxplot()+
    facet_grid(.~id,scales = "free_y")
#### End ####

#### Calprotectin vs activity ####
# meta%>%filter(cal_group=="Yes"|cal_group=="No",!is.na(activity_category))%>%
#     ggplot(aes(x=activity_category,y=cal_value))+
#     geom_point()
  
drug %>%
   group_by(library_id) %>%
   summarize(drug_cate = list(drug_cate))%>%
    filter(!is.na(drug_cate))%>%
  ggplot(aes(x = drug_cate)) +
       geom_bar() +
       scale_x_upset()

# 
# 
# g<-tmp%>%as_tibble()%>%left_join(drug)%>%
# ggplot(aes(x = drug_cate,fill=infection)) +
#   geom_bar() +
#   scale_x_upset()
  #   
  # pdf("output/main.figure/anelloviridae.drug_group.pdf", width = 10, height = 7)
  # g
  # dev.off()
  
df<-dna.filter%>%filter(family=="Anelloviridae")
df<-df%>%inner_join(meta)
#df$infection<-"Uninfected"
#df[df$cov>=0.33,"infection"]<-"Infected"
tmp<-df %>% filter(immu_suppress=="Yes",study_group=="VEO") 

g<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Anti-TNF")%>%
  ggplot(aes(x = activity_category,y=read)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point()
pdf("output/main.figure/anti-tnf.severity.rpkm.pdf", width = 7, height = 5)
g
dev.off()

test<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Anti-TNF")%>%select(read,activity_category)
  kruskal.test(read ~ activity_category,data=test)
  summary(aov(read ~ activity_category,data=test))
test<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Steroids")%>%select(read,activity_category)
  kruskal.test(read ~ activity_category,data=test)
  summary(aov(read ~ activity_category,data=test))
test<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Inhibitor")%>%select(read,activity_category)
  kruskal.test(read ~ activity_category,data=test)
  summary(aov(read ~ activity_category,data=test))


g<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Steroids")%>%
  ggplot(aes(x = activity_category,y=read)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point()
pdf("output/main.figure/inhibitor.severity.rpkm.pdf", width = 7, height = 5)
g
dev.off()

g<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Inhibitor")%>%
  ggplot(aes(x = activity_category,y=read)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point()
pdf("output/main.figure/steroids.severity.rpkm.pdf", width = 7, height = 5)
g
dev.off()


g1<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Anti-TNF")%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar(position = "fill") 
g2<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Anti-TNF")%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar() 

pdf("output/main.figure/anti-tnf.severity.virus.pdf", width = 7, height = 5)
g1
g2
dev.off()

g1<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Steroids")%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar(position = "fill") 
g2<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Steroids")%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar() 

pdf("output/main.figure/steroids.severity.virus.pdf", width = 7, height = 5)
g1
g2
dev.off()

g1<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Inhibitor")%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar(position = "fill") 
g2<-tmp%>%as_tibble()%>%left_join(drug)%>%filter(drug_cate=="Inhibitor")%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar() 

pdf("output/main.figure/inhibitor.severity.virus.pdf", width = 7, height = 5)
g1
g2
dev.off()



tmp<-df %>% filter(!is.na(activity_category),study_group=="VEO") %>% select(library_id,activity_category,infection)

g<-tmp%>%
  ggplot(aes(x = activity_category,fill=infection)) +
  geom_bar() 

pdf("output/main.figure/anelloviridae.drug_group.activity.pdf", width = 10, height = 7)
g
dev.off()

#### End ####

#### Extra ####
df<-dna.filter%>%filter(family=="Anelloviridae")
df<-df%>%inner_join(meta)
df$infection<-"Uninfected"
df[df$cov>=0.33,"infection"]<-"Infected"
tmp<-df %>% filter(study_group=="VEO"|study_group=="Healthy") %>% dplyr::select(study_group,infection)
tmp$study_group<-droplevels(tmp$study_group)
fisher.test(table(tmp))[1][[1]]


tmp<-df %>% filter(immu_suppress=="No",cal_group=="Yes"|cal_group=="No") %>% select(cal_group,infection)
tmp$cal_group<-droplevels(tmp$cal_group)
fisher.test(table(tmp))[1][[1]]

#### End ####

#### Ecoli O157 ####
read_reads<-function(x){
  p = read.delim(x, header=F,
                 stringsAsFactors = FALSE)
  p<-p%>%mutate(rpkm=V3/(V2/1000*(sum(V3)+sum(V4)/1000000)))
  names(p)<-c("Genome", 
              "Genome.len",
              "mapped", 
              "unmapped",strsplit(strsplit(x,"\\/")[[1]][length(strsplit(x,"\\/")[[1]])],"\\.")[[1]][1])
  
  return(p[-nrow(p),c(1,5)])
}

path = "/media/lorax/users/guanxiang/2_veo/liang.veo/input/decontam"
filename <- list.files(path = path,full.names = TRUE,pattern = ".ecoli$") 
l<-c()

l<-lapply(filename,read_reads)
refseq<-Reduce(function(x, y) left_join(x, y,by="Genome"),l)

tmp<-refseq%>%
  gather(key=library_id,value=ecoli,-Genome) %>%
  left_join(meta)

ggplot(tmp%>%filter(study_group=="VEO"|study_group=="Healthy"),aes(x=study_group,y=ecoli*100,color=study_group))+
  geom_boxplot()


wilcox.test((tmp%>%filter(study_group=="VEO"))$ecoli,(tmp%>%filter(study_group=="Healthy"))$ecoli)
#### End ####












