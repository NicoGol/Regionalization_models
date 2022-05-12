ncovr <- st_read('NAT.shp')

if (TRUE) {
  library(sf)
  library(rgeoda)
  guerry_path <- system.file("extdata", "Guerry.shp", package = "rgeoda")
  guerry <- st_read(guerry_path)
  queen_w <- queen_weights(guerry)
  data <- guerry[c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids')]
  guerry_clusters <- redcap(4, queen_w, data, "fullorder-completelinkage")
  guerry_clusters
}
