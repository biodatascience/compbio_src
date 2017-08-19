library(rvest)
library(readr)
orgs <- read_tsv("organisms.tsv")
res <- lapply(orgs$species, getdetails)
details <- do.call(rbind, res)
details
stopifnot(orgs$species == details$org)
orgs <- cbind(orgs, details[,-1])
write_tsv(orgs, "organisms_genes.tsv")

getdetails <- function(org) {
  readit <- function(org, prefix) {
    read_html(paste0("http://",prefix,".ensembl.org/",
                     org, "/Info/Annotation"))
  }
  page <- tryCatch(readit(org, "www"),
                   error=function(e) readit(org, "plants"))
  tabs <- page %>% html_nodes("table") %>% html_table()
  cutafterspace <- function(x) sub("(.*?) .*","\\1",x)
  getit <- function(tab,x) {
    text <- tab[grep(x,tab[,1])[1],2]
    as.numeric(gsub(",","",cutafterspace(text)))
  }
  offset <- if (colnames(tabs[[1]])[2] == "Number of gene models") 1 else 0
  bp <- getit(tabs[[1+offset]], "Golden Path Length")
  coding <- getit(tabs[[2+offset]], "Coding genes")
  noncoding <- getit(tabs[[2+offset]], "Non coding genes")
  pseudo <- getit(tabs[[2+offset]], "Pseudogenes")
  tx <- getit(tabs[[2+offset]], "Gene transcripts")
  data.frame(org, bp, coding, noncoding, pseudo, tx)
}
