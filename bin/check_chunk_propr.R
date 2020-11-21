
blocks2combs <- function(counts, nblocks){
  
  ngroup <- rep(1:nblocks, each = 100)
  leftover <- ncol(counts) %% nblocks
  if(leftover > 0) ngroup <- c(ngroup, rep(nblocks + 1, leftover))
  split <- split(1:ncol(counts), ngroup)
  
  combs <- expand.grid(1:length(split), 1:length(split))
  combs <- t(apply(combs, 1, sort))
  combs <- unique(combs) 
  combs <- combs[combs[,1] != combs[,2],]
  
  return(list(combs, split))
}

propr.check <- function(counts, rho){

  nblocks = ncol(counts) %/% 100

  l <- blocks2combs(counts, nblocks)
  combs <- l[[1]]
  split <- l[[2]]

  diffs = c()
  for(i in 1:nrow(combs)){
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]

    rho.i = propr(counts, "rho", p=20, select=c(batch1, batch2))
    rho.j = subset(rho, select=c(batch1, batch2))
    if(identical(rho.i@matrix, rho.j@matrix)==FALSE){
      diffs = c(diffs, i)
    }
  }
  print(diffs)
}





