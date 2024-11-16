# Package
version     = "0.6.1"
author      = "Charles Blake"
description = "Self-contained fit of linear models with regression diagnostics"
license     = "MIT/ISC"

# Deps
requires    "nim >= 1.6.0"
requires    "cligen >= 1.7.8"
requires    "spfun >= 0.7.3"
skipDirs    = @["fitl"]
installExt  = @[".nim"]
bin         = @["fitl", "fitl/qtl", "fitl/gof", "fitl/dists", "fitl/estMI"]
