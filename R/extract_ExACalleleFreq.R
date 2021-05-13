extract_ExACalleleFreq <- function(rowData){


  # using a POST call perform a bulk query of ExAC database using the variant names
  res<- httr::POST("http://exac.hms.harvard.edu/rest/bulk/variant/variant", body = jsonlite::toJSON(rowData$ExACname))

  # convert the data from JSON format to list
  dat<-httr::content(res, "text") %>% jsonlite::fromJSON()

  dat %>% head()
}
