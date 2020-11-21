
# library paths
.libPaths(c("/nfs/software/R/packages","/nfs/users2/cn/ierb/R/", "/nfs/users2/cn/sjin/R/x86_64-redhat-linux-gnu-library/3.3/"))


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
  make_option(c("-t", "--tissue"),
              type="character",
              default=NULL,
              help="A file containing the tissue name, or just the tissue name.",
              metavar="character"),
  make_option(c("-p", "--prefix"),
              type="character",
              default=NULL,
              help="Output prefix. Eg. prefix=tissue3, then output: tissue3_results.csv, tissue3_fdr.csv, etc.",
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
              help="Interval for updateCutoff step")
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

# read tissue
if (file.exists(opt$tissue))
  {
    tissue = readLines(opt$tissue)  # eg. "Brain - Cerebellar Hemisphere"
    if(is.null(opt$prefix)){
      basetissue = basename(opt$tissue)  # eg. tissue2
    }else{
      basetissue = opt$prefix   
    }
  }else{
    tissue = opt$tissue
    basetissue = opt$prefix
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


# ======= #
# TISSUES #
# ======= #

# get tissue samples
tis=which(colData(rse_gene)[,"smtsd"]==tissue)
tissamp=colData(rse_gene)[tis,"run"]

# get the tissue samples that have a reasonable expression
Mt=t(m[mygenes,tissamp])

# reduce dataset, if test
if (opt$test) {
  set.seed(0)
  smp = sample(1:ngenes, opt$test_size, replace=FALSE)
  Mt = Mt[, smp]
}

# zero replacement
M=cmultRepl(Mt,method="CZM",label=0,output="p-counts")


# ====== #
# MEMORY #
# ====== #

# free variables to save memory
rm(m)
rm(av)
rm(mygenes)
rm(tis)
rm(tissamp)
rm(Mt)
rm(rse_gene)  


# =============== #
# PROPORTIONALITY #
# =============== #

# run proportionality analysis
pro=propr(M,metric="rho", p=20)
pro=updateCutoffs(pro,seq(0.4,0.95,opt$cutoff_interval))

# output 
if (!dir.exists(opt$outdir)) {dir.create(opt$outdir, recursive=TRUE)}

# write FDR
out1 = file.path(opt$outdir, paste(basetissue, "_fdr.csv", sep=""))
write.table(pro@fdr, out1, row.names = FALSE, quote=FALSE, sep=",", dec=".")

# write pair results
out2 = file.path(opt$outdir, paste(basetissue, "_results.csv", sep=""))
write.table(pro@results, out2, row.names = FALSE, quote=FALSE, sep=",", dec=".")