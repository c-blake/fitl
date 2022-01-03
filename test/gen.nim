import std/[math, random, strformat], spfun/[gauss, cauchy]

proc gen(x0=0.0, dx=0.25, n=40, p=2, mix=0.0, gaussScale=1.0, cauchScale=1.0,
         reseed=false) =
  ## Generate data to test linear model fitting robustness.  Polynomial is:
  ## `y = sum j=0..p x^j/(j+1) + Gaussian & Cauchy noise mixture`.
  if reseed:
    randomize()
  var x = x0
  for i in 0 ..< n:
    var y = if rand(1.0) >= mix: gaussScale*gauss.qtl(rand(1.0))
            else               : cauchScale*cauchy.qtl(rand(1.0))
    for j in 0..p: y += float(1)/float(j+1)*x^j
    stdout.write &"{y:.8g}"
    for j in 0..p: stdout.write &" {x^j:.7g}"
    stdout.write '\n'
    x += dx

when isMainModule:
  import cligen
  dispatch gen, help={
    "x0"        : "starting x",
    "dx"        : "x step",
    "n"         : "num points",
    "p"         : "polynomial order",
    "mix"       : "fraction of Cauchy samples",
    "gaussScale": "scale of Gaussian noise",
    "cauchScale": "scale of Cauchy noise",
    "reseed"    : "randomize RNG seed"}
