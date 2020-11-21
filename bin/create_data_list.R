# ======== #
# SETTINGS #
# ======== #

# load library path
.libPaths(c("/nfs/software/R/packages","/nfs/users2/cn/ierb/R/", "R/x86_64-redhat-linux-gnu-library/3.3/"))

# load library
library(recount)


# ========= #
# LOAD DATA #
# ========= #

# parameters
args = commandArgs(trailingOnly=TRUE)
data = args[1]  # input data file
out.folder = args[2]

# read data
if (file.exists(data))
  {
    load(file.path(data))
  }else{
    download_study(data, type = "rse-gene")
    load(file.path("/users/cn/sjin/projects/proportionality/data", data, "rse_gene.Rdata"))
  }


# ======= #
# TISSUES #
# ======= #

# get all tissues
tissues=unique(colData(rse_gene)[,"smtsd"])

# get number of samples in each tissue
nsam1=rep(0,length(tissues))
for (i in 1:length(tissues)){
	nsam1[i]=length(colData(rse_gene)[which(colData(rse_gene)[,"smtsd"]==tissues[i]),"sampid"])
}

# tissues with >= 20 samples: 50
mytissues=tissues[which(nsam1>20)]


# ======== #
# PATIENTS #
# ======== #

# get all donors
donors=unique(substring(colData(rse_gene)[,"sampid"],1,9))

# get number of samples from each donor 
nsam2=rep(0,length(donors))
for (i in 1:length(donors)){
	nsam2[i]=length(colData(rse_gene)[which(substring(colData(rse_gene)[,"sampid"],1,9)==donors[i]),"sampid"])
}

# donors with >= 20 samples
# (note that some donors might have various replicate samples from the same tissue)
mydonors=donors[which(nsam2>=20)]


# =========== #
# OUTPUT LIST #
# =========== #

# create output folders
for (i in c("donors", "tissues"))
  {
    o2 = file.path(out.folder, i)
    if (!dir.exists(o2 )) {dir.create(o2)}
  }

# output donors
i = 1
for (donor in mydonors)
  {
    write( donor, file.path(out.folder, "donors", paste("donor", i, sep="")) )
    i = i + 1
  }

# output tissues
i = 1
for (tissue in mytissues)
  {
    write( tissue, file.path(out.folder, "tissues", paste("tissue", i, sep="")) )
    i = i + 1
  }