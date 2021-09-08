mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("HMDB00812","HMDB00425","HMDB03609","HMDB00191","HMDB00300","HMDB00211","HMDB00133","HMDB00225","HMDB00694","HMDB03290","HMDB00195","HMDB01241","HMDB00086","HMDB60256","HMDB00848","HMDB01385","HMDB01564","HMDB02092")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_2_", "net", "png", 72, width=NA)
