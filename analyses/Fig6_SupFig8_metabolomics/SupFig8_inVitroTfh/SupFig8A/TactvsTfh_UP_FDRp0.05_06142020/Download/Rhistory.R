mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("HMDB00190","HMDB06469","HMDB12880","HMDB00162","HMDB00187","HMDB05175","HMDB02250","HMDB00210","HMDB00624","HMDB03976","HMDB00123","HMDB00208","METPA0491","HMDB00262","HMDB06471","HMDB01565","HMDB00043","HMDB05066","HMDB01553","HMDB00097","HMDB06317","HMDB01173","HMDB00019","HMDB00273","HMDB00126","HMDB06210","HMDB00791","HMDB00099","HMDB02366","HMDB00254","HMDB00824","HMDB00156","HMDB02013","HMDB00491","HMDB00828","HMDB13336","HMDB00134","HMDB00062","HMDB13127","HMDB00052","HMDB00201")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "hmdb");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 72, width=NA)
