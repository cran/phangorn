## zzz.R 

.packageName <- "phangorn"


.onLoad  <- function(libname, pkgname) {
    library.dynam("phangorn", pkgname, libname)
}

