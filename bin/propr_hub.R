# ========= #
# ARGUMENTS #
# ========= #

library(optparse)

option_list = list(
  make_option(c("-i", "--input"), 
              type="character", 
              default=NULL, 
              help="Input file containing proportionality pairs.", 
              metavar="character"),
  make_option(c("-c", "--cutoff"),
              type="character",
              default=NULL,
              help="File containing cutoff information."),
  make_option(c("-o", "--output"),
              type="character",
              default=NULL,
              help="Output file with hubs stats",
              metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# ========= #
# LIBRARIES #
# ========= #

library(igraph)


# ======== #
# ANALYSIS #
# ======== # 

print("----reading and processing input file")

# read input 
csv = read.csv(opt$input)
head(csv)

# cutoff
cutoff = round(min(csv[,"propr"]), 2)
print(cutoff)

# proportionality pairs
ppairs = csv[,1:2]
nppairs = nrow(ppairs)

print("----getting graph")

# graph
graph=graph_from_data_frame(ppairs,directed=FALSE)

print("----getting graph degrees")

# graph degree
deg=degree(graph)
ndeg = sum(deg)
ncomp = length(deg)

print("----getting hubs")

# hubs with >= 10, or >=100
hub10 = deg[deg>=10]
nhub10 = length(hub10)
hub100 = deg[deg>=100]
nhub100 = length(hub100)
hub1000 = deg[deg>=1000]
nhub1000 = length(hub1000)
hub10000 = deg[deg>=10000]
nhub10000 = length(hub10000)

print("----print output")

d = data.frame('cutoff'=cutoff, 'ppairs'= nppairs, 'ncomp' = ncomp,'ndeg'= ndeg, 
                'min'= min(deg), 'q25'= as.numeric(summary(deg)[2]), 
                'median'= median(deg), 'mean'= round(mean(deg),2), 
                'q75'= as.numeric(summary(deg)[5]), 'max'= max(deg),
                'hub10'= nhub10, 'hub100'= nhub100, 'hub1000'=nhub1000, 'hub10K'=nhub10000)
print(d)

write.table(d, opt$output, row.names = FALSE, quote=FALSE, sep=",", dec=".")

# print(paste("CUTOFF: ", cutoff, sep=""))
# print(paste("PPAIRS: ", nppairs, sep=""))
# print(paste("DEG: ", ndeg, sep=""))
# print(paste("MIN: ", min(deg), sep=""))
# print(paste("Q25: ", as.numeric(summary(deg)[2]), sep=""))
# print(paste("MEDIAN: ", median(deg), sep=""))
# print(paste("MEAN: ", mean(deg), sep=""))
# print(paste("Q75: ", as.numeric(summary(deg)[4]), sep=""))
# print(paste("MAX: ", max(deg), sep=""))
# print(paste("HUB10: ", nhub10, sep=""))
# print(paste("HUB100: ", nhub100, sep=""))