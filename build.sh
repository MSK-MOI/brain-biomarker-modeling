R CMD INSTALL bbutilities
R CMD build bbutilities
R CMD check bbutilities_1.0.0.tar.gz
R CMD Rd2pdf bbutilities