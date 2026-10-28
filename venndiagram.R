
library(VennDiagram)

vennlist <- get_vennlist(obj=COVID2_agggenus, factorNames="Group")
vennp <- venn.diagram(vennlist,
height=5,
width=5,
filename=NULL,
fill=c("red", "green", "blue"),
)
grid::grid.draw(vennp)
