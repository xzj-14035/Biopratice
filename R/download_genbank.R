#' @title  download_genbank
#' @description This function downloads genbank files from NCBI.
#' @param acc The accession number
#' @param db The database, it can be nuccore, nucest, nucgss or popset
#' @param destfile The path where you want to save the genbank files
#' @export
download_genbank <- function(acc, db, destfile){
  db <- db
  db <- paste('db=', db, '&id=', sep = '')
  accs <- acc
  base_url <- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
  for (acc in accs){
    mature_url <- paste(base_url, db, acc, '&rettype=gb', sep = '')
    gb.file_name <- paste(acc, '.gb', sep = '')
    destfile <- paste(destfile, gb.file_name, sep = '/')
    utils::download.file(mature_url,destfile = destfile)
  }
}
