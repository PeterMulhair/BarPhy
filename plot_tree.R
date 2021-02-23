#!/usr/bin/env Rscript

library(ggtree)
library(getopt)

#Get call input
spec <- matrix(c('phylo', 't', 1, "character", 'output', 'o', 1, "character"), byrow=TRUE, ncol=4)

opt <- getopt(spec)

#Read rooted tree file
tree <- read.tree(opt$phylo)

#Create tree object
if(tree$Nnode<80) {
p <- ggtree(tree) + geom_tiplab(offset = 0.00001, aes(color=grepl("query", label)))
} else {p <- ggtree(tree, layout="circular") + geom_tiplab2(aes(color=grepl("query", label)))}


#Write tree to pdf
if(tree$Nnode<40){
pdf(opt$output, "_tree.pdf", height = 20, width = 50)
} else { if(tree$Nnode<60) {
pdf(opt$output, "_tree.pdf", height = 50, width = 50)
} else {pdf(opt$output, "_tree.pdf", height = 150, width = 150)
}
}
print(p)
dev.off()
