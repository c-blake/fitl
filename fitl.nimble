# Package
version     = "0.2.0"
author      = "Charles Blake"
description = "Self-contained fit of linear models with regression diagnostics"
license     = "MIT/ISC"

# Deps
requires    "nim >= 1.6.0"
requires    "cligen >= 1.5.24"
requires    "spfun >= 0.3.0"
installExt  = @[ ".nim" ]
bin         = @[ "fitl", "fitl/qtl" ]
