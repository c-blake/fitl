# Package
version     = "0.6.4"
author      = "Charles Blake"
description = "Self-contained fit of linear models with regression diagnostics"
license     = "MIT/ISC"

# Deps
requires    "nim >= 1.6.0"
requires    "cligen >= 1.9.2"
requires    "spfun >= 0.7.6"
skipDirs    = @["fitl"]
installExt  = @[".nim"]
bin         = @["fitl", "fitl/qtl", "fitl/gof", "fitl/dists", "fitl/estMI"]
