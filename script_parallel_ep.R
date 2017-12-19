
---
  title: ""
author: ""
date: '22 ноября 2017 г '
output: pdf_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r clearence, R libraries}
rm(list = ls())

library(ape)
library(reldna)
```


Reading in E. coli K12 genome (GenBank Accesion U00096.2), creating reverse complement genome and tranformation genome into character form 
```{r Reading in E. coli genome}
e.coli_U00096.2_char<-read.GenBank(access.nb = 'U00096.2', as.character = T)$U00096.2
e.coli_U00096.2_string <- paste0(e.coli_U00096.2_char, collapse = '')[[1]]



shifted_by <- 200
extended_e.coli_U00096.2_string <- paste0(
  substr(e.coli_U00096.2_string, nchar(e.coli_U00096.2_string)-shifted_by, nchar(e.coli_U00096.2_string)),
  e.coli_U00096.2_string,
  substr(e.coli_U00096.2_string, 0, shifted_by)
)

```

```{r Parallel}
extended_ecoli_string <- extended_e.coli_U00096.2_string
#pseudo_tss <- 1:250000

library(parallel)
library(reldna)
no_cores = detectCores() - 2
cl <- makeCluster(no_cores)

calculate_EP_on_interval <- function(tss_position, extended_string, shifted_by, zout) {
  p <- lseqspline1D(substr(extended_string, tss_position, tss_position+shifted_by+150), bound=c(50, 350), ref=251 )
  return (p$mpot[p$x %in% zout])
}

clusterExport(cl, 'extended_ecoli_string')
clusterExport(cl, 'calculate_EP_on_interval')
clusterEvalQ(cl, {library(reldna)})

#res <- parSapply(cl, X = pseudo_tss, FUN = function(x) calculate_EP_on_interval(x, extended_ecoli_string, 250, -480:239))
#stopCluster(cl)

#save(res, file = '/home/jane/Документы/Misha/mol_a_2018/parsapply_ep_ecoli1strand_1_.rda')
```


```{r By fragments}
bounds <- c(seq(1, nchar(extended_ecoli_string), 250000), nchar(extended_ecoli_string))

for (i in 1:(length(bounds) - 1)){
  pseudo_tss <- bounds[i]:bounds[i+1]
  print(length(pseudo_tss))
  
  res <- parSapply(cl, X = pseudo_tss, FUN = function(x) calculate_EP_on_interval(x, extended_ecoli_string, 250, -480:239))
  
 # assign( paste0('sliced_ep_', min(pseudo_tss), '_', max(pseudo_tss)), res)
  
  rm(list = 'res')
  save( res, file = paste0('/home/jane/Документы/Misha/mol_a_2018/',  paste0('sliced_ep_', min(pseudo_tss), '_', max(pseudo_tss)), '.rda'))
}


```


