items=unique(ra_items$column)
plot_list=list()
for (a in items){
x=pheatmap(as.matrix(data),
             cluster_rows = F,
             fontsize = 14,
             main = a)
  plot_list[[a]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}
g<-do.call(grid.arrange,plot_list)
ggsave("g.pdf",g)