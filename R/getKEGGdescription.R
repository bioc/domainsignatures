## Map KEGG Ids to description of pathways. Essentially, this is a wrapper
## around the KEGGPATHID2NAME environment.

getKEGGdescription <- function(ids)
  data.frame(ID=ids, description=unlist(mget(substring(ids, 4,100),
                       KEGGPATHID2NAME)))

