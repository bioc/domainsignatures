## ==========================================================================
## The show method for 'ipDataSource' objects
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("show", signature="ipDataSource",
          definition=function(object){
            cat("object of class 'ipDataSource' containing mapping data",
                "\nfor type '", object@type, "' pathways\n", sep="")
            cat("   ", length(object@pathways), " pathways with " ,
                length(object@genes), " genes\n   and ", length(object@domains),
                " annotated InterPro domains\n", sep="")
          }
          )
