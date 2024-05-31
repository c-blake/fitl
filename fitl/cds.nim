when not declared(stdin): import std/[syncio, formatfloat]
import std/[strutils, strformat, algorithm, random]
from fitl/qtl    import quantile
from spfun/binom import initBinomP, est

proc cds*(x: seq[float], m=50, sort=false): seq[seq[float]] =
  ## Use an *already sorted* `x` to make `m` Parzen Qmid re-samples, maybe sort.
  result.setLen m
  for i in 0..<m:
    result[i].setLen x.len
    for j in 0..<x.len:
      result[i][j] = x.quantile(rand(1.0))
    if sort:
      result[i].sort

proc cdswarm*(iput="", oput="rs", m=20, gplot="", ci=0.95) =
  ## Read numbers from stdin|`iput` if not-""; Emit Parzen-interpolated quantile
  ## samples (of the *same* sample size) to files `{oput}NNN`.
  var x: seq[float]
  for f in lines(if iput.len>0: iput.open else: stdin): x.add f.strip.parseFloat
  x.sort; let n = x.len
  let xs = x.cds(m, sort=true)
  let g = if gplot.len > 0: open(gplot, fmWrite) else: nil
  let e = open(&"{oput}EDF", fmWrite)
  let tagL = &"{0.5 - 0.5*ci:.03f}"
  let tagH = &"{0.5 + 0.5*ci:.03f}"
  let l = open(&"{oput}{tagL}" , fmWrite)
  let h = open(&"{oput}{tagH}" , fmWrite)
  for j, f in x:
    e.write f, "\n"
    let (lo, hi) = initBinomP(j, n).est(ci)     #XXX Check alignment & maybe
    l.write f," ",lo,"\n"                       #..  emit leading & trailing
    h.write f," ",hi,"\n"                       #..  edges to connect to 0,1
  e.close
  if g != nil: g.write &"""#set terminal png size 1920,1080 font "Helvetica,10"
#set output "rsSwarm.png"
set key top left noautotitle    # EDFs go bot left->up right;Dot keys crowd plot
set style data steps
set xlabel "Sample Value"
set ylabel "Probability"
set linetype 1 lc rgb "blue" lw 3
set linetype 2 lc rgb "red"  lw 1
set linetype 3 lc rgb "red"  lw 2
set linetype 4 lc rgb "black" dashtype 0
plot """
  for i, x in xs:
    let opath = &"{oput}{i:03}"
    let o = open(opath, fmWrite)
    for f in x: o.write f, "\n"
    o.close
    if g != nil:
      g.write (if i==0: "" else: ",\\\n     "), &"'{opath}' u 1:($0/{n}) ls 4"
  if g != nil:
    g.write &",\\\n     '{oput}EDF' u 1:($0/{n}) title 'EDF' ls 1"
    g.write &",\\\n     '{oput}{tagL}' title 'EDF{tagL}' ls 2"
    g.write &",\\\n     '{oput}{tagH}' title 'EDF{tagH}' ls 3"
    g.write "\n"; g.close

when isMainModule:
  when not declared(stdin): import std/[syncio, formatfloat]
  when defined danger: randomize()
  import cligen; dispatch cdswarm, help={
    "iput" : "input path or \"\" for stdin",
    "oput" : "output path prefix; outs Get numbered",
    "m"    : "number of resamples; e.g. for plots",
    "gplot": "generate a gnuplot script to plot",
    "ci"   : "CI for Wilson score confidence bands"}
