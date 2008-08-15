## ===========================================================================
## ipDataSource
## ---------------------------------------------------------------------------
## A class to store mapping information between genes, pathways and
## interPro domains
## ---------------------------------------------------------------------------

setClass("ipDataSource",
         representation(genes="character", pathways="character",
                        domains="character", gene2Domains="list",
                        path2Domains="list", dims="numeric",
                        type="character"))

