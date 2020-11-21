
read_data <- function(data) {
    if (file.exists(data))
    {
        load(file.path(data))
    }else{
        download_study(data, type = "rse-gene")
        load(file.path("/users/cn/sjin/projects/proportionality/data", data, "rse_gene.Rdata"))
    }
}

