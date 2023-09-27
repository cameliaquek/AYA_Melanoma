#######################################
# Sankey and drug groups
#######################################

# Load library
library(tidyverse) 
library(networkD3) 
library(ggplot2)


# Set up data frame
dat <- read.delim(file = "~/links_setup.txt")

links <- data.frame( 
  source=dat[,1], 
  target=dat[,2], 
  value=dat[,3] ) 

nodes <- data.frame( 
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique() 
) 
nodes %>% head(3) 

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1 
links %>% head(3) 


# Draw and save Sankey plot
p <- sankeyNetwork(Links = links, Nodes = nodes, height=800,width=800,iterations = 100,
              colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20);"),
              Source = "IDsource", Target = "IDtarget",Value = "value",
              NodeID = "name",
              #LinkGroup = "LinkGroupID",
              #NodeGroup = "NodeGroupID",
              fontSize= 16, 
              nodeWidth=30,nodePadding=10,margin=-0.12,
              sinksRight=F)
              
pdf(file = "treatment-switch_sankey.pdf")
p 
dev.off()

