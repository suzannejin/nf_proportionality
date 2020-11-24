# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option(c("-i", "--input"), 
              type="character", 
              default=NULL, 
              help="Input dataset file. If only the ID is given, then the corresponding dataset will be downloaded.", 
              metavar="character"),
  make_option(c("-o", "--outdir"),
              type="character",
              default=getwd(),
              help="Output directory. Default = current working directory [%default].",
              metavar="character"),
  make_option(c("-d", "--donor"),
              type="character",
              default=NULL,
              help="A file containing the donor name, or just the donor name.",
              metavar="character"),
  make_option(c("-p", "--prefix"),
              type="character",
              default=NULL,
              help="Output prefix. Eg. prefix=donor3, then output: donor3_results.csv, donor3_fdr.csv, etc.",
              metavar="character"),
  make_option(c("-l", "--datadir"),
              type="character",
              default="/users/cn/sjin/projects/proportionality/data",
              help="Directory where the dataset will be downloaded, if not provided by user. Default=%default",
              metavar="character"),
  make_option("--test",
              action="store_true",
              default=FALSE,
              help="Only test on a small subset of data."),
  make_option("--test_size",
              type="integer",
              default=200,
              help="Size of the test data. Default = %default",
              metavar="number"),
  make_option("--ncores",
              type="integer",
              default=1,
              help="Number of cpus to compute the updateCutoffs step.",
              metavar="number"),
  make_option("--cutoff_interval",
              type="double",
              default=0.01,
              help="Interval for updateCutoff step"),
  make_option("--interval_min",
              type="double",
              default=0.3,
              help="Minimum cutoff for FDR computation"),
  make_option("--interval_max",
              type="double",
              default=0.95,
              help="Maximum cutoff for FDR computation"),
  make_option("--chunk_size",
              type="integer",
              default=100,
              help="Size for propr.chunk",
              metavar="number"),
  make_option("--local",
              action="store_true",
              default=FALSE,
              help="Load local libraries.")
); 
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

if (opt$local){
  # library paths
  .libPaths(c("/nfs/software/R/packages","/nfs/users2/cn/ierb/R/", "/nfs/users2/cn/sjin/R/x86_64-redhat-linux-gnu-library/3.3/"))
}

library(recount)
library(zCompositions)
library(propr)
library(doMC)


# ========= #
# FUNCTIONS #
# ========= #

blocks2combs <- function(counts, n){

  # define blocks
  nblocks = ncol(counts) %/% n
  ngroup <- rep(1:nblocks, each = n)
  leftover <- ncol(counts) - length(ngroup)
  if(leftover > 0) ngroup <- c(ngroup, rep(nblocks + 1, leftover))

  # check size
  if (length(unique(ngroup)) <= 2){
    stop(paste("ERROR: chunk size ", n, " is too big for data frame of size [", nrow(counts), "][", ncol(counts), "]", sep=""))
  }

  # split groups
  split <- split(1:ncol(counts), ngroup)
  combs <- expand.grid(1:length(split), 1:length(split))
  combs <- t(apply(combs, 1, sort))
  combs <- unique(combs) 
  combs <- combs[combs[,1] != combs[,2],]

  print(paste("----runing propr for ", nrow(combs), " chunks of size ", n))

  return(list(combs, split))
}

chunk2propr <- function(i, chunk, metric="rho", ivar=NA, symmetrize=FALSE,
                        alpha, p=100, interval=seq(0.4,0.95,0.01) ){

  print(i)

  # compute propr chunk
  rho.i <- propr(chunk, metric = metric, ivar = ivar, alpha = alpha, p=p)
  rho.i <- updateCutoffs(rho.i,interval)

  # get cutoff
  cutoff <- min(rho.i@fdr[,"cutoff"][which(rho.i@fdr[,"FDR"]<.05)])

  # update interval and recompute, if no reasonable cutoff
  icut = 1
  while(cutoff==Inf && icut<5){
    m <- max(interval)
    b <- (interval[2]-interval[1])/10
    print(paste("Re-updating cutoff using a different interval (", m, ",", 1-b, ",", b, ")", sep=""))
    interval <- seq(m, 1-b, b)
    rho.i <- updateCutoffs(rho.i,interval)
    cutoff <- min(rho.i@fdr[,"cutoff"][which(rho.i@fdr[,"FDR"]<.05)])
    icut = icut+1
  }
  if(cutoff==Inf){
    m <- min(rho.i@fdr[,"FDR"])
    stop(paste("No cutoff for FDR < 0.05, minimum = ", m, sep=""))
  }

  # get fdr at cutoff
  fdr <- rho.i@fdr[,"FDR"][which(rho.i@fdr[,"cutoff"]==cutoff)]
  print(paste("Obtained cutoff[", cutoff, "] and fdr[", fdr, "] for chunk ", i, sep=""))

  return(list(rho.i@matrix, cutoff, fdr))
}

chunk2full <- function(counts, RES, split, combs){
  
  # define variables
  QUILT <- matrix(0, ncol(counts), ncol(counts))
  cutoff <- 0.0
  fdr <- 0.0

  # collect chunks
  for(i in 1:nrow(combs)){

    # update cutoff and fdr
    cutoff <- cutoff + RES[[i]][[2]]
    fdr <- fdr + RES[[i]][[3]]

    # Fill final matrix with each job
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]
    patch.i <- c(batch1, batch2)
    QUILT[patch.i, patch.i] <- RES[[i]][[1]]
  }

  # average cutoff and fdr
  cutoff <- cutoff / i
  fdr <- fdr / i

  # rename columns & rows
  matrix <- QUILT
  rownames(matrix) <- colnames(counts)
  colnames(matrix) <- colnames(counts)

  return(list(matrix, cutoff, fdr))
}

file2res <- function(RES, combs, dir){

  print(paste("----collecting ", nrow(combs), " chunk files previously saved in ", dir, sep=""))

  for(i in 1:nrow(combs)){
    file2 <- paste0(dir, "/job-", combs[i,1], "+", combs[i,2], ".csv")
    csv <- read.csv(file2, row.names = 1, header= TRUE)
    matrix <- data.matrix(csv)
    RES[[i]][[1]] <- matrix
  }

  return(RES)
}

propr.chunk <- function(counts, metric = c("rho", "phi", "phs", "cor", "vlr"),
                        ivar = "clr", symmetrize = FALSE, alpha, p=100,
                        n=100, ncores = 1, interval=seq(0.3,0.95,0.01),
                        dir = NA){
  
  # divide data into chunks
  l <- blocks2combs(counts, n)
  combs <- l[[1]]
  split <- l[[2]]
  
  # parallelize propr computation for i chunks
  doMC::registerDoMC(cores = ncores)
  `%dopar%` <- foreach::`%dopar%`
  
  RES <- foreach::foreach(i = 1:nrow(combs)) %dopar% {

    # get chunk
    batch1 <- split[[combs[i,1]]]
    batch2 <- split[[combs[i,2]]]
    chunk = subset(clr, select = c(batch1, batch2))

    # compute propr 
    l_propr <- chunk2propr(i, chunk, metric=metric, ivar=ivar, symmetrize=symmetrize, alpha=alpha, p=p, interval=interval)

    if(is.na(dir)){
      # return propr matrix, cutoff, and concrete FDR
      l_propr
    }else{
      # save data if required
      file2 <- paste0(dir, "/job-", combs[i,1], "+", combs[i,2], ".csv")
      print(paste("----saving tmp file[", i, "] to ", file2, sep=""))
      write.csv(l_propr[[1]], file=file2)

      # return cutoff and FDR
      list(NULL, l_propr[[2]], l_propr[[3]])
    }
  }

  # collect files
  if(!is.na(dir)){
    RES <- file2res(RES, combs, dir)
  }

  # merge chunks
  l_propr <- chunk2full(counts, RES, split, combs)

  return(l_propr)
}


# ========= #
# LOAD DATA #
# ========= #

print("----loading data")

# read donor
if (file.exists(opt$donor))
  {
    donor = readLines(opt$donor)  # eg. GTEX-OHPN
    if(is.null(opt$prefix)){
      basedonor = basename(opt$donor)
    }else{
      basedonor = opt$prefix   # eg. donor3
    }
  }else{
    donor = opt$donor
    basedonor = opt$prefix
  }

# read gtex
if (file.exists(opt$input))
  {
    load(file.path(opt$input))
  }else{
    download_study(opt$input, type = "rse-gene")
    load(file.path(opt$datadir, opt$input, "rse_gene.Rdata"))
  }


# ============ #
# FILTER GENES #
# ============ # 

print("----filtering genes")

#data matrix
## gene names as rows, and sample names as columns
## each column-row combination has one integer value
m=assay(rse_gene)

#simlify the gene names (the ones with dots aren't recognized)
rownames(m)=substring(rownames(m),1,15)

#average rows (genes) across columns (samples)
## so we get the average expression by gene
av=apply(m,1,mean)

#leave out the lower 25% of genes for a reasonable expression
mygenes=rownames(m)[which(av>boxplot(av)$stats[2,1])]  
ngenes = length(mygenes)

#we should replace zeros right here, but it's too slow: 
#M=cmultRepl(t(m[mygenes,]),method="CZM",label=0,output="p-counts")


# ====== #
# DONORS #
# ====== #

print("----getting donor matrix")

# get samples from the donor
don=which(substring(colData(rse_gene)[,"sampid"],1,9)==donor)
donsamp=colData(rse_gene)[don,"run"]

# get the donor samples that have a reasonable expression
Md=t(m[mygenes,donsamp])

# reduce dataset, if test
if (opt$test) {
  set.seed(0)
  smp = sample(1:ngenes, opt$test_size, replace=FALSE)
  Md = Md[, smp]
}

# replace zeros
M=cmultRepl(Md,method="CZM",label=0,output="p-counts")


# ====== #
# MEMORY #
# ====== #

# free variables to save memory
rm(m)
rm(av)
rm(mygenes)
rm(don)
rm(donsamp)
rm(Md)
rm(rse_gene)


# =============== #
# PROPORTIONALITY #
# =============== #

print("----running proportionality analysis")

# output folder
if (!dir.exists(opt$outdir)) {dir.create(opt$outdir, recursive=TRUE)}
dir2 = paste(opt$outdir, "tmp", sep="/")
if (!dir.exists(dir2)) {dir.create(dir2, recursive=TRUE)}

# CLR transformation
clr = propr:::proprCLR(M)

#ch = propr.chunk(clr, "rho", ivar=NA, alpha=NA, ncores=3, p=20, n=400, interval=seq(.3, .95, .01), dir=dir2)

# run proportionality analysis per chunk 
ch = propr.chunk(clr, metric="rho", ivar=NA, alpha=NA, ncores=opt$ncores, p=20, n=opt$chunk_size,
                 interval=seq(opt$interval_min, opt$interval_max, opt$cutoff_interval), dir=dir2)
matrix = ch[[1]]
cutoff = ch[[2]]
fdr = ch[[3]]

# organize results data frame
labels <- propr:::labRcpp(ncol(clr))
# lrv <- propr:::lr2vlr(clr)   # somehow this step does not work for me
results <-
    data.frame(
      "Partner" = labels[[1]],
      "Pair" = labels[[2]],
      # "lrv" = propr:::lltRcpp(lrv),
      # "metric" = factor(metric),
      # "alpha" = factor(alpha),
      "propr" = propr:::lltRcpp(matrix)
    )

# get proportional pairs
mypairs=which(results[,"propr"]>=cutoff)
results = results[mypairs,]


# ====== #
# OUTPUT #
# ====== #

# proportionality pairs
out1 = file.path(opt$outdir, paste(basedonor, "_results.csv", sep=""))
write.table(results, out1, row.names = FALSE, quote=FALSE, sep=",", dec=".")

# cutoff
out2 = file.path(opt$outdir, paste(basedonor, "_cutoff.txt", sep=""))
writeLines(as.character(cutoff), out2)

# fdr
out3 = file.path(opt$outdir, paste(basedonor, "_fdr.txt", sep=""))
writeLines(as.character(fdr), out3)

# remove tmp files
unlink(dir2, recursive=TRUE)