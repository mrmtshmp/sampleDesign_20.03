#' Make data
#' PI Dr.Endoh
#' 2020/09/14

dir.sub <- './sub'
Bibtex <- TRUE

list.files.dir.sub <- list.files(path = dir.sub)

for(i in 1:length(list.files.dir.sub))
  source(sprintf("%s/%s", dir.sub, list.files.dir.sub[i]))


colinfo <- 
  readxl::read_excel(
    path = sprintf("%s/%s",dir.data,fn.data),
    sheet = 2
    )

data <- 
  readxl::read_excel(
    path = sprintf("%s/%s",dir.data,fn.data),
    col_names = colinfo$col_names,
    col_types = colinfo$col_types,
    skip=1
    ) %>%
  data.frame()

save(data,file = sprintf("%s/%s", dir.ADS, fn.ADS))


