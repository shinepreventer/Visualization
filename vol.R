rm(list = ls())
library(ggrepel)
library(ropls) 
library(ggplot2) 
data = read.table("met.txt", header=T, row.names=1, sep="\t")
group = read.table("group.txt", header=T, row.names=1, sep="\t")
plsda = opls(t(data), group$group, predI = 2,orthoI=0)
vol0 = data.frame(plsda@vipVn,stringsAsFactors = FALSE)
colnames(vol0) = "Vip"
A.mean = apply(data[,1:6],1,FUN = mean)
B.mean = apply(data[,7:10],1,FUN = mean)
FC = B.mean/A.mean
vol=vol0
vol$FC = FC
log2FC = log(FC,2)
vol$log2FC = log2FC
pvalue = apply(data, 1, function(x)
{t.test(x[1:6],x[7:10])$p.value})
p.log = -log10(pvalue)
vol$logp = p.log
vol$p = pvalue
vol$type <- 'insig'
vol$type[vol$FC>=2 & vol$p <= 0.05 & vol$Vip>1] <- 'up'
vol$type[vol$FC<=0.5 & vol$p <= 0.05 & vol$Vip>1] <- 'down'
vol$A = "no"
vol$A[vol$type =='up'| vol$type =='down'] = 'sig'
vol$comp = rownames(vol)
write.csv(vol,"vol.csv",row.names = F)
ggplot(vol, aes(log2FC,logp)) 
  geom_point(aes(fill = type,size=Vip),shape = 21 )+
  scale_fill_manual(values = c("#66c2a5", "grey", "#f46d43"))+
  geom_vline(xintercept = c(-1, 1), lty=3, lwd=0.8,color = "black") + 
  geom_hline(yintercept = -log10(0.05), lty=3,color = 'red', lwd=0.8) +
  geom_text_repel(aes(label = comp), data = subset(vol, type %in% c('up', 'down')), max.overlaps = 15) + 
  labs(title="",
       x = 'Log2FC',
       y = '-log10 pvalue')+
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.ticks.length = unit(0.15, "cm"), 
        axis.line = element_line(color = "black"), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        legend.text = element_text(color = 'black', size = 12)) 
  
ggsave("vol.pdf",units = "in",width = 15,height = 8,dpi = 1200)

