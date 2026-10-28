#isntallation
devtools::install_github("david-barnett/microViz@0.9.3")

#run library
library(microViz)

#Interactive ordination exploration (pseq is your phyloseq object)
ord_explore(pseq)