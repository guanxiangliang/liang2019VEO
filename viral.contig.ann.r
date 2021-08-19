uniprottax<-read.delim("/media/lorax/users/guanxiang/1_igram/liang2019/input/database/uniprot.virome/uniprot.virome.meta",header=F)
names(uniprottax)<-c("uniprot.hit.id","taxid","name")
uniprottax$taxid<-as.numeric(as.character(uniprottax$taxid))
pa="/media/lorax/users/guanxiang/2_veo/liang.veo/input/cross.assembly/"
## RNA contig 
rna.contig<-rna.contig[rowSums(rna.contig[,80:84]<25)>=4,] # remove contamination contigs
rna.ann<-rna.contig%>%select(contig)
datapath <- file.path(paste(pa,"rna",".orf.number",sep=""))
orfn<-read.delim(datapath,sep="\t",header=F)
orfn<-separate(data = orfn, col = V1, into = c("contig"), sep = "\\_")
tmp<-data.frame(table(orfn$contig))
names(tmp)<-c("contig","orf.number")
rna.ann<-inner_join(rna.ann,tmp)
# Uniprot data process
datapath <- file.path(paste(pa,"rna",".uniprot",sep=""))
uniprot<-read.delim(datapath,header=F,sep="\t")
#uniprot<-read.delim("./veo.uniprot.2",header=F,sep="\t")
names(uniprot)<-c("contig", "uniprot.hit.id", "pident", "qlen","slen","alignment.len" ,"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
uniprot<- uniprot  %>% 
  separate(col = contig, into = c("contig","orfs.id"), sep = "\\_") 
uniprot<-inner_join(uniprot,uniprottax,by="uniprot.hit.id")
uniprot<-cbind(uniprot,getTaxonomy(uniprot$taxid,"/media/lorax/users/guanxiang/1_igram/liang2019/input/database/ac2tax/accessionTaxa.sql"))
# Uniprot orf number
tmp<-data.frame(table(uniprot$contig))
names(tmp)<-c("contig","uniprot.total.orf.number")
rna.ann<-inner_join(rna.ann,tmp,by="contig")
# Calculate uniprot species
tmp<-uniprot[,c("contig","species")]
tmp<-table(tmp)
tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
tmp<-uniprot[,c("contig","species")]
tmp<-data.frame(table(tmp))
tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
tmp3<-cbind(tmp1,tmp2)
names(tmp3)<-c("uniprot.species","uniprot.species.orf.number")
tmp3$contig<-rownames(tmp3)
rownames(tmp3)<-NULL
rna.ann<-inner_join(rna.ann,tmp3,by="contig")
rna.ann$uniprot.species<-as.character(rna.ann$uniprot.species)
rna.ann[rna.ann$uniprot.species.orf.number==0, ][, "uniprot.species"] <- "Others"
# Calculate uniprot genus
tmp<-uniprot[,c("contig","genus")]
tmp<-table(tmp)
tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
tmp<-uniprot[,c("contig","genus")]
tmp<-data.frame(table(tmp))
tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
tmp3<-cbind(tmp1,tmp2)
names(tmp3)<-c("uniprot.genus","uniprot.genus.orf.number")
tmp3$contig<-rownames(tmp3)
rownames(tmp3)<-NULL
rna.ann<-inner_join(rna.ann,tmp3,by="contig")
rna.ann$uniprot.genus<-as.character(rna.ann$uniprot.genus)
rna.ann[rna.ann$uniprot.genus.orf.number==0, ][, "uniprot.genus"] <- "Others"
# Calculate uniprot family
tmp<-uniprot[,c("contig","family")]
tmp<-table(tmp)
tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
tmp<-uniprot[,c("contig","family")]
tmp<-data.frame(table(tmp))
tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
tmp3<-cbind(tmp1,tmp2)
names(tmp3)<-c("uniprot.family","uniprot.family.orf.number")
tmp3$contig<-rownames(tmp3)
rownames(tmp3)<-NULL
rna.ann<-inner_join(rna.ann,tmp3,by="contig")
rna.ann$uniprot.family<-as.character(rna.ann$uniprot.family)
rna.ann[rna.ann$uniprot.family.orf.number==0, ][, "uniprot.family"] <- "Others"

# Calculate uniprot Order
tmp<-uniprot[,c("contig","order")]
tmp<-table(tmp)
tmp1<-data.frame(colnames(tmp)[apply(tmp,1,which.max)])
tmp<-uniprot[,c("contig","order")]
tmp<-data.frame(table(tmp))
tmp2<-data.frame(tapply(tmp$Freq, tmp$contig, max))
tmp3<-cbind(tmp1,tmp2)
names(tmp3)<-c("uniprot.order","uniprot.order.orf.number")
tmp3$contig<-rownames(tmp3)
rownames(tmp3)<-NULL
rna.ann<-inner_join(rna.ann,tmp3,by="contig")
rna.ann$uniprot.order<-as.character(rna.ann$uniprot.order)
rna.ann[rna.ann$uniprot.order.orf.number==0, ][, "uniprot.order"] <- "Others"
rna.ann<-merge(rna.ann,host[,c("family","host.name")],by.x="uniprot.family",by.y="family",all.x=T)
rna.ann<-rna.ann%>%filter(grepl("RNA",host.name))
rna.ann<-rna.ann%>%mutate(ratio=uniprot.total.orf.number/orf.number)%>%filter(ratio>0.5)