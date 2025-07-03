library(randomForest)
library(rfPermute) 
library(caret)
otu = read.csv("2.csv",row.names = 1,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE) 
dim(otu)
head(otu)

spe = otu
spe[3:ncol(spe)] <- sweep(spe[3:ncol(spe)],1,rowSums(spe[3:ncol(spe)]),'/')*100
spe
classification <- factor(spe[[2]])

categories <- unique(classification)
num_categories <- length(categories)

#spe$grazing = as.factor(spe$grazing)
set.seed(12345)
RF = randomForest(x = spe[-c(1:2)],
                  y = factor(spe[[2]]), 
                  mtry= floor(sqrt(ncol(spe[-c(1:2)]))),
                  ntree = 500, 
                  importance=TRUE,
                  localImp = TRUE,	
                  proximity = TRUE,
                  oob.prox = TRUE,
                  norm.votes = TRUE,
                  na.action=na.omit, 
                  keep.inbag=TRUE,
                  sampsize= nrow(spe[-c(1:2)]),
                  nodesize = 1,
                  maxnodes = NULL, 
                  nPerm=1,
                  do.trace=FALSE, 
                  corr.bias=TRUE,
                  xtest=NULL,
                  ytest = NULL,
                  weights=NULL, 
                  replace = TRUE, 
                  classwt = NULL, 
                  cutoff = rep(1 / num_categories, num_categories), 
                  strata=NULL,
)

print(RF)
#summary(RF)
RF$call 
RF$type 
RF$classes 
RF$y 
RF$ntree 
RF$mtry 
RF$inbag 

RF$predicted
RF$forest 

set.seed(12345);
rf = randomForest(x = spe[-c(1:2)],
                  y= factor(spe[[2]]),
                  importance=TRUE,
                  ntree = 10000
)

#rf$err.rate 
which.min(rf$err.rate[,1])
rf$err.rate[which.min(rf$err.rate[,1]),] 
#summary(rf) 
plot(rf) 
legend("bottomright",legend = colnames(rf$err.rate),pch=20,col = c(1:4))

detect_ntree = function(i,...){
  set.seed(12345);
  rf = randomForest(x=spe[-c(1:2)],y=factor(spe[,2]),
                    importance=TRUE,
                    ntree = i,...);
  res = data.frame(ntree = i,class.error = mean(rf$confusion[,4]));
  return(res);
}

ntree = sapply(seq(1,500,1), 
               detect_ntree
)
#ntree
ntree[,which.min(ntree[2,])]

maxnode = data.frame(row.names =seq(3,nrow(spe),1))
for (i in seq(3,nrow(spe),1)){
  set.seed(12345);
  rf = randomForest(x = spe[-c(1:2)],
                    y= factor(spe[[2]]),
                    maxnodes = i,
                    importance=TRUE,
                    ntree = 415
  )
  maxnode[paste(i,"",sep = ""),1] = mean(rf$confusion[,4])
}

library(dplyr)
maxnode = arrange(maxnode,V1) 
maxnode 


detect_mtry = function(i){
  set.seed(12345);
  rf = randomForest(x = spe[-c(1:2)], y= factor(spe[[2]]),
                    maxnode=14,
                    importance=TRUE,
                    mtry =i,ntree = 415)
  res = data.frame(mtry = i,class.error = mean(rf$confusion[,4]))
  #res = data.frame(mtry = i,error.rate = rf$err.rate[nrow(rf$err.rate),1])
  return(res)
}

mtry = sapply(seq(10,ncol(spe[-c(1:2)]),10), 
              detect_mtry
)
#mtry
mtry[,which.min(mtry[2,])]


names = c()
for (i in name){
  names = c(names,paste(i,seq(2,50,1),sep = "_"))
}
length(names) 
err = lapply(strsplit(names,"_"),as.numeric) 
detect_para = function(err){
  #library(randomForest);
  set.seed(12345); 
  rf = randomForest(x = spe[-c(1:2)],
                    y= factor(spe[[2]]),
                    importance=TRUE,
                    maxnodes = err[[1]],
                    mtry = err[[2]],
                    ntree = err[[3]]                  
  )
  res = data.frame(maxnode = err[[1]],
                   mtry = err[[2]],
                   ntree = err[[3]], # class error
                   class.error = mean(rf$confusion[,4]) ,
                   ntree = which.min(rf$err.rate[,1]),
                   error.rate = rf$err.rate[which.min(rf$err.rate[,1]),1]
  ) 
  return(res)
}

library(parallel);
cores = detectCores();
cl = makeCluster(cores);
clusterExport(cl, c("spe"));
clusterEvalQ(cl, library(randomForest));
system.time(
  {
    results <- parSapply(cl, err, detect_para); 
  } 
)
stopCluster(cl);
results[,which.min(results[4,])] 



otu = read.csv("2.csv",row.names = 1,header = TRUE, quote = "", comment.char="",check.names = FALSE,stringsAsFactors = F) # 微生物组数据
dim(otu)
head(otu)


spe = otu
spe[3:ncol(spe)] <- sweep(spe[3:ncol(spe)],1,rowSums(spe[3:ncol(spe)]),'/')*100
head(spe)

library(randomForest)
set.seed(12345) 
RF.best = randomForest(x = spe[-c(1:2)], 
                       y= factor(spe[,2]),
                       importance=TRUE,
                       maxnodes = 14,
                       mtry = 160,
                       ntree =415,
                       Perm = TRUE,
                       )   #class.error=0.1428571

rfPermute::confusionMatrix(RF.best) 
summary(RF.best) 
imp = data.frame(importance(RF.best),MDA.p = RF.best$importanceSD[4])
head(imp)

library(dplyr)
imp = arrange(imp,desc(MeanDecreaseGini)) 
head(imp) 
write.csv(imp,"importance.csv",quote = FALSE)
#install.packages("tibble")
library(tibble)
library(dplyr)

compound_ids <- read.csv("compound_ids.csv", header = TRUE, stringsAsFactors = FALSE)
head(compound_ids)
compound_ids <- compound_ids %>% dplyr::rename(CompoundID = CompoundID)
impc <- imp %>% rownames_to_column(var = "CompoundID")
filtered_imp <- inner_join(compound_ids, impc, by = "CompoundID")
head(filtered_imp)
rownames(filtered_imp) <- filtered_imp$CompoundID
filtered_imp <- filtered_imp[, -which(names(filtered_imp) == "CompoundID")]
head(filtered_imp)
filtered_imp = arrange(filtered_imp,desc(MeanDecreaseGini)) 
write.csv(filtered_imp, "filtered_importance_with_rownames.csv", quote = FALSE)
## top10
imp10 = filtered_imp[1:10,]
imp10
library(ggplot2)
p1= ggplot(imp10,aes(x=MeanDecreaseGini,y=reorder(rownames(imp10),MeanDecreaseGini)))+
  geom_bar(position = position_dodge(),
           width = 0.5,
           stat = "identity",
           fill="steelblue")+ 
  theme_minimal() +
  xlab("Mean Decrease in Gini Index")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(size = 16,colour = "black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.x.bottom = element_text(size=16,color="black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #panel.grid = element_blank()
  )
p1

library(rstatix)
library(dplyr)
all= read.csv("all.csv",row.names = 1,header = TRUE, quote = "", comment.char="",check.names = FALSE,stringsAsFactors = F) 
dim(all)
head(all)
all[3:ncol(all)] <- sweep(all[3:ncol(all)],1,rowSums(all[3:ncol(all)]),'/')*100
BC = all[c("depth",rownames(imp10))]
BC

sig = matrix(data=NA,nrow = 10,ncol = 2,dimnames = list(names(BC[2:11]),c("species","sig")))
for (i in 2:11){
  assign(paste(names(BC[i]),"kw",sep="."),kruskal.test(BC[,i],as.factor(BC[,1])));
  tmp = get(paste(names(BC[i]),"kw",sep="."))
  sig[i-1,] = c(names(BC[i]),ifelse(tmp$p.value >0.05,"ns"," *"));
} # assign()将kruskal.test()
kw = mget(paste(names(BC[2:11]),"kw",sep=".")) 
capture.output(kw,file = "kw.list.txt",append = FALSE)
sig = as.data.frame(sig)
head(sig)
BC1 = data.frame(name = rownames(BC),BC[c(1:11)])
BC1 
BC1 = gather(BC1,key="species","abundance",-name,-depth)
BC1
BC2 = BC1 %>% group_by(depth,species) %>% 
  get_summary_stats(abundance,type = "mean_se") 
BC2

da = merge(BC2,sig,by.x = c("species"),by.y = c("species"))
da$group = factor(da$depth,levels = c("A","B","C"))
da[da$depth %in% c("A","B"),"sig"] = NA 
da2 = merge(da,imp,by.x = c("species"),by.y = "row.names")
da2

p2 = ggplot(da2,
            aes(x=mean,y=reorder(species,MeanDecreaseGini),fill=depth))+
  geom_bar(position = position_dodge(),
           width = 0.5,
           stat = "identity")+
  geom_errorbar(aes(xmin=mean-se,xmax=mean+se),
                width=0.3,
                color="black",linewidth=0.5,
                position = position_dodge(0.5))+
  geom_text(aes(x=(mean+se),label=sig,group=depth),
            hjust =0,nudge_x=0.5, 
            nudge_y=0.075, 
            angle = -90, 
            size=4.75,fontface="bold")+
  theme_minimal() +
  xlab("Abundance")+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(axis.text.y = element_blank(), 
        #axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=14,color="black"),
        axis.title.x.bottom = element_text(size=16,color="black"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank() 
  )
p2 

library(patchwork)

ppp <- p1+p2+plot_layout(nrow = 1, widths = c(2, 2),
                  guides = "collect")
ppp
