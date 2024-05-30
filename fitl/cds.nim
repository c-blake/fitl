when not declared(stdin): import std/[syncio, formatfloat]
import std/[strutils, strformat, algorithm, random]
from fitl/qtl import quantile

proc cds*(x: seq[float], m=50, sort=false): seq[seq[float]] =
  ## Use an *already sorted* `x` to make `m` Parzen Qmid re-samples, maybe sort.
  result.setLen m
  for i in 0..<m:
    result[i].setLen x.len
    for j in 0..<x.len:
      result[i][j] = x.quantile(rand(1.0))
    if sort:
      result[i].sort

proc cdswarm*(iput="", oput="rs", m=50, sort=true, gplot="") =
  ## Read numbers from stdin|`iput` if non-empty; Emit Parzen-interpolated
  ## quantile samples (of the *same* sample size) to files `{oput}NNN`.
  var x: seq[float]
  for f in lines(if iput.len>0: iput.open else: stdin): x.add f.strip.parseFloat
  x.sort
  let xs = x.cds(m, sort)
  let gp = if gplot.len > 0: open(gplot, fmWrite) else: nil
  if gp != nil: gp.write """#set terminal png size 1920,1080 font "Helvetica,10"
#set output "rsHerd.png"
set style data steps; set key noautotitle       # key/legend crowds out plot
set linetype 1 lc rgb "black"; set linetype 2 lc rgb "black"
set linetype 3 lc rgb "black"; set linetype 4 lc rgb "black"
set linetype 5 lc rgb "black"; set linetype 6 lc rgb "black"
set linetype 7 lc rgb "black"; set linetype 8 lc rgb "black"
set linetype cycle 8                            # cycle *must be* >= 8
plot"""
  for i, x in xs:
    let opath = &"{oput}{i:03}"
    let o = open(opath, fmWrite)
    for f in x: o.write f, "\n"
    o.close
    if gp != nil: gp.write (if i==0: " " else: " ,"), &"'{opath}' u 1:0"
  if gp != nil: gp.close

when isMainModule:
  when not declared(stdin): import std/[syncio, formatfloat]
  when defined danger: randomize()
  import cligen; dispatch cdswarm, help={
    "iput" : "input path or \"\" for stdin",
    "oput" : "output path prefix; outs Get numbered",
    "m"    : "number of resamples; e.g. for plots",
    "sort" : "sort the resamples to easy plotting",
    "gplot": "generate a gnuplot script to plot"}
