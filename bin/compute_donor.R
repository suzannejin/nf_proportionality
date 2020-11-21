
# library paths
# .libPaths(c("/nfs/software/R/packages","/nfs/users2/cn/ierb/R/", "/nfs/users2/cn/sjin/R/x86_64-redhat-linux-gnu-library/3.3/"))


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
  make_option("--only_results",
              action="store_true",
              default=FALSE,
              help="Only compute proportionality, without calculating FDR.")
); 
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ========= #
# LIBRARIES #
# ========= #

library(recount)
library(zCompositions)
library(propr)


# ========= #
# LOAD DATA #
# ========= #

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

# run proportionality analysis
pro=propr(M, metric="rho", p=20)

# output 
if (!dir.exists(opt$outdir)) {dir.create(opt$outdir, recursive=TRUE)}

# write pair proportionality
out1 = file.path(opt$outdir, paste(basedonor, "_results.csv", sep=""))
write.table(pro@results, out1, row.names = FALSE, quote=FALSE, sep=",", dec=".")

if (opt$only_results==FALSE){
  
  # Calculate FDR
  pro=updateCutoffs(pro,seq(0.4,0.95,opt$cutoff_interval))  # ncore > 1 don't work

  # write FDR
  out2 = file.path(opt$outdir, paste(basedonor, "_fdr.csv", sep=""))
  write.table(pro@fdr, out2, row.names = FALSE, quote=FALSE, sep=",", dec=".")

}



