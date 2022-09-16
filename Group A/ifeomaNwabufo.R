#merging tables
ex=read.csv("data/exKIRC.csv")
ds=read.csv("data/dsKIRC.csv")
ds=subset(ds,select=c(Drug.name,Cell.line.name,IC50)) #selecting the non-null columns
ex<-t(ex)#transposing the dataframe
colnames(ex)<-ex[1,] #reassigning the column names
ex<-ex[-1,]
rownames(ex) <- sub(pattern = ".", "-", row.names(ex), fixed = TRUE)#Converting "." to "-" in row names of ex data frame 
ex <- cbind('Cell.line.name' = rownames(ex), ex)
df = merge(x = ds, y = ex, by = "Cell.line.name")

#Mdensity plots code
library(ggExtra)
library(tidyverse)
library(ggplot2)
data=read.csv("merged(1).csv")
geneIndex = 1
drug = "CGP-60474"
df = subset(data, select = c(Cell.line.name,Drug.name,IC50, A1BG ) )%>% filter(Drug.name=="CGP-60474")
med=median(df$A1BG)
df$GeneExpressLevel = ifelse (df$A1BG >= med, "high", "low")
x=ggplot(df, aes(x=A1BG, y=IC50, label=Cell.line.name))+
  geom_point(size=2, aes(colour=GeneExpressLevel))+
  geom_point(size=2)+
  geom_smooth(method=('lm'))+
  scale_colour_manual(values=c('darkorange', 'grey54'))+
  geom_text(nudge_x=0, nudge_y=0.2, size=6, colour='darkcyan')+
  theme_bw()+
  theme(text=element_text(size=20), legend.position='bottom')
ggMarginal(x,type="density",margins="both",groupColour=TRUE)