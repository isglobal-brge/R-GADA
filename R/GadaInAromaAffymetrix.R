addGadaToAromaAffymetrix <-  function()
  {
    #Make sure aroma.affymetrix is installed...

    require('aroma.affymetrix') || throw("Package not loaded: aroma.affymetrix");
    
    source(system.file("exec/addGadaToAromaAffymetrix.R", package = "gada"));

  }
