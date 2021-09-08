mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("C00379","C00183","C00299","C01620","C00199","C00474","C00134","C00802","C00407","C00294","C00262","C00387","C00623","C00025","C02341","C00165","C05422","C00327","C00099","C00049","C06104","C00020","C00212")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "kegg");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 72, width=NA)
