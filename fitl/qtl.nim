proc parzen*[F:SomeFloat, U:SomeFloat](x: openArray[F], q: U): U =
  ##[ Return Parzen Qmid Quantile of sorted openArray.  Ma, Genton &Parzen 2011:
  "Asymptotic properties of sample quantiles of discrete distributions" has more
  discussion.  Personally, I see this as better than any alternative definition
  & also the right generalization of ancient "mid-ranking ties" ideas in rank
  correlation calculations (or even Wilcoxon's original 1945 paper).  Qmid is
  sadly not as widely used|even known as it should be. ]##
  let n = x.len; let N = float(n)
  if n < 1: return 0.0
  let qN = float(q) * N
  if qN <= 0.5: return x[0]
  if qN >= N - 0.5: return x[n - 1]
  var iL, jL, iH, jH: int
  var xL, xH, cH, cL: float
  jL = int(qN); iL = jL                   #XXX Could replace scans with exponen.
  xL = x[iL]                              #  ..expanding search for *many* dups.
  while iL > 0 and x[iL-1] == xL: dec(iL)
  while jL < n and x[jL] == xL: inc(jL)
  cL = 0.5 * float(iL + jL)
  if cL <= qN:                            #Scan higher
    iH = jL; jH = iH
    xH = if jH < n: x[jH] else: xL
    while jH < n and x[jH] == xH: inc(jH)
    cH = 0.5 * float(iH + jH)
  else:                                   #Scan lower
    cH = cL; iH = iL; jH = jL; xH = xL    # *H <- *L
    jL = iH; iL = jL
    xL = if iL > 0: x[iL - 1] else: xH
    while iL > 0 and x[iL - 1] == xL: dec(iL)
    cL = 0.5 * float(iL + jL)
  let r = (qN - cL) / (cH - cL)
  return U((1.0 - r) * xL  +  r * xH)

import spfun/beta, fitl/basicLA
proc harrellDavis*[F:SomeFloat](w: var seq[F], n: int, q: F, err=1e-14) =
  ##[ Set weights `w` on `n` order stats to estimate quantile `q`. Harrell-Davis
  1982: A New Distribution-Free Quantile Estimator has details.  `err` is the
  convergence error for betaI func eval which should be small because weights
  are differences of said across successive index vals.  For the same reason,
  internal calc is in `float`, though weights may be stored as `float32`. ]##
  w.setLen n
  let (p, q) = (float(n + 1)*q, float(n + 1)*(1 - q))
  let nInv = 1.0/float(n)
  var last = betaI(0.0, p, q, err)
  for i in 0..<n:
    let curr = betaI(float(i + 1)*nInv, p, q, err)
    w[i] = curr - last
    last = curr

type Method* = enum Parzen, HarrellDavis

proc quantile*[F:SomeFloat, U:SomeFloat](x: openArray[F], q: U, m=Parzen,
                                         w: ptr F=nil): U =
  case m
  of Parzen: parzen x, q
  of HarrellDavis:
    if w.isNil:
      var w: seq[F]; harrellDavis w, x.len, q
      dot x[0].addr, w[0].addr, x.len
    else:
      dot x[0].addr, w, x.len

when isMainModule:
  when not declared(stdin): import std/[syncio, formatfloat]
  import std/[strutils, algorithm, random], cligen
  when defined danger: randomize()
  proc qtl(`method`=Parzen, ps: seq[float]) =
    ## Read one column of numbers on stdin; Emit Parzen-interpolated quantiles
    ## for probabilities `ps`. Eg.: *echo 1 1 2 4|tr ' ' '\\n'|qtl 0 .5 1*
    ## writes *1.0 1.66.. 4.0*.
    ##
    ## If `ps[0]<0` & only `p` (eg. ``-210``) instead emit -THIS many re-samples
    ## to stdout.  Eg. `split(1)` can then turn into 10 files for 21 data.
    var x: seq[float]
    for line in stdin.lines: x.add parseFloat(line.strip)
    x.sort
    if ps.len == 1 and ps[0] < 0:
      for i in 1 .. -int(ps[0]): echo x.quantile(rand(1.0))
    else:
      for i, p in ps:
        if i > 0: stdout.write " "
        stdout.write x.quantile(p, `method`)
      stdout.write '\n'
  include cligen/mergeCfgEnv;dispatch qtl,help={"method":"Parzen, HarrellDavis"}
