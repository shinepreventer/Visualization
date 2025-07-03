library(tidyverse)
library(phyloseq)
library(EasyStat)
library(ggClusterNet)
library(tidyverse)
library(ggClusterNet)
library(phyloseq)
library(igraph)
library(ggalt)
library(ggnewscale)
library(gtable)
library(grid)
library(igraph)
library(sna)
library(tidyfst)
library(grid)
library(gtable)
library(igraph)
library(network)
library(SpiecEasi)
library(WGCNA)


metadata = read.table("./metadata.tsv", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
otutab = read.table("./otu_table.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
taxonomy = read.table("./taxonomy.tsv", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

# Extract only those ID in common between the two tables
idx = rownames(otutab) %in% rownames(taxonomy)
otutab = otutab[idx,]
taxonomy = taxonomy[rownames(otutab),]


ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy)))
ps


netpath = "all" 
dir.create(netpath) 
set.seed(123)



all =  network.pip(ps = ps,
                   N=0,
                   big = FALSE,
                   select_layout = FALSE,
                   #layout_net = "model_maptree2",
                   r.threshold=0.70,
                   p.threshold=0.05,
                   group = "all",
                   label = FALSE,
                  # fill = "Phylum",
                   #size = "igraph.degree",
                   #path = netpath ,
                   zipi = FALSE,
                   method = "sparcc",
                   method = "pearson",
                   ncpus = 80
)

netpath= "./resultall"
dir.create(netpath)

allnode=allpip[[2]][["net.cor.matrix"]][["node"]]
alledge=allpip[[2]][["net.cor.matrix"]][["edge"]]
write.csv(allnode,paste0(path,"/allnode.csv"))
write.csv(alledge,paste0(path,"/alledge.csv"))

saveRDS(all,"cor.matrix.all.group.rds")



