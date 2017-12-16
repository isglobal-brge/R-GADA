
############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("gada", pkg, lib)
}

.onAttach <- function (lib, pkg){
  data(genomicInfo)  ##RPR automatic load of the genome info...  
}

.onUnload <- function(libpath)
    library.dynam.unload("gada", libpath)


############ End of .First.lib ###############


