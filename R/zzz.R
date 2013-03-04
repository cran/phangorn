## zzz.R 

.packageName <- "phangorn"


.onLoad  <- function(libname, pkgname) {
    library.dynam("phangorn", pkgname, libname)
}

.aamodels <- c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU","Blossum62","Dayhoff_DCMut","JTT-DCMut")
