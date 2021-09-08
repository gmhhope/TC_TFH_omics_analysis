mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c("C00379","C00105","C01620","C00474","C00013","C00134","C00346","C00149","C00130","C00294","C00258","C00122","C02341","C05422","C06104","C00020","C00212","C00055")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "kegg");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetSMPDB.PathLib(mSet, "hsa")
mSet<-SetOrganism(mSet, "hsa")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
mSet<-PlotPathSummary(mSet, "path_view_0_", "png", 72, width=NA)
mSet<-PlotPathSummary(mSet, "path_view_0_", "svg", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
