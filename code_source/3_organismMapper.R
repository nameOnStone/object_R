organismMapper <- function(x){
  if(x=="human"){
    orgdb <- "org.Hs.eg.db"
    species <- "hsa" 
  } else if (x=="mouse"){
    orgdb <- "org.Mm.eg.db"
    species <- "mmu" 
  } else if(x=="rat"){
    orgdb <- "org.Rn.eg.db"
    species <- "rno" 
  }
  return(list("orgdb"=orgdb,"species"=species))
}