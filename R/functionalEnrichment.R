.enrichmentTest = function(targets,background,targetsFunc,backgroundFunc)
{
  ft = fisher.test(cbind(c(targetsFunc,backgroundFunc),c(targets,background)),alternative = "greater")
  return(c(ft$estimate,ft$conf.int[1],ft$conf.int[2],ft$p.value,targetsFunc,backgroundFunc))
}

.enrichmentSet = function(targets,background,set)
{
  genesInSet = names(set)[set==1]
  backgroundFunc = intersect(background,genesInSet)
  targetsFunc = intersect(backgroundFunc,targets)
  c(.enrichmentTest(length(targets),length(background)-length(targets),length(targetsFunc),length(backgroundFunc)-length(targetsFunc)),length(genesInSet))
}

.enrichmentTable = function(targets,background,table)
{
  out = matrix(data = 0,nrow=ncol(table),ncol=8,dimnames = list(colnames(table),c("Odds ratio","conf_int_min","conf_int_max","p.value","bonf p.value","No. targets in set","No. targets in background","Total set size")))
  for(s in 1:ncol(table))
  {
    out[s,c(1:4,6:8)] = .enrichmentSet(targets,background,table[,s])
  }
  out[,5] = p.adjust(out[,4],method = "bonferroni")
  out[order(out[,5]),]
}

#' @export
functionalEnrichment = function(targets,background)
{
  modernOut = .enrichmentTable(targets,background,table=modern)
  wormcat1 = .enrichmentTable(targets,background,table=wormcat1)
  wormcat2 = .enrichmentTable(targets,background,table=wormcat2)
  wormcat3 = .enrichmentTable(targets,background,table=wormcat3)
  return(list("modern" = modernOut,"WormCat1" = wormcat1,"WormCat2"=wormcat2,"WormCat3"=wormcat3))
}
