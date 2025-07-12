const adx = defined(release) and defined(useAdix)
when adx: import adix/nsort             # Usually faster linear-time number sort
else    : import std/algorithm
proc rsort*[E](xs: var seq[E]) = ## ~ radix `sort` & `sorted`
  when adx:
    when E is float32: nsort xs, 0u32, xfFlt
    else:              nsort xs, 0u64, xfFlt # !float32 => float64
  else: sort xs
proc rsorted*[E](xs: openArray[E]): seq[E] = result = xs; rsort result

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

proc walshAverages*[F: SomeFloat](x: openArray[F]): seq[F] =
  ##[ Walsh 1950: Note on a theorem due to Hoeffding; Ann.Math.Stat. 21(1) first
  framed discussion of pairwise-including-self averages. ]##
  let n = x.len; var k = 0
  result.setLen n*(n + 1)div 2
  for i in 0..<n:                         # Self-inclusive pairing important
    for j in 0..i: result[k] = 0.5*(x[i] + x[j]); inc k

proc hodgesLehmann*[F:SomeFloat, U:SomeFloat](x: openArray[F], q: U): U =
  ##[ Return Generalized Hodges-Lehmann Estimator Parzen-quantile(pairwise
  averages) of sorted openArray `x`, but extreme quantiles not recommended. ]##
  var pairAvs = walshAverages x
  rsort pairAvs
  parzen pairAvs, q

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

proc hDW*[F:SomeFloat, U:SomeFloat](x: openArray[F], q: U): U =
  ##[ Return q-Generalized Hodges-Lehmann Estimator (pairwise averages) of
  sorted openArray `x`, but with Harrell-Davis rather than Parzen quantiles. ]##
  var pAs = walshAverages x
  rsort pAs
  var w: seq[F]; harrellDavis w, pAs.len, q
  dot pAs[0].addr, w[0].addr, pAs.len

type Method* = enum Parzen="parzen", HarrellDavis="d", HodgesLehmann="l",
                    HDW="b" # b for both

proc quantile*[F:SomeFloat, U:SomeFloat](x: openArray[F], q: U, m=Parzen,
                                         cbα=0.05, w: ptr F=nil): U =
  case m
  of Parzen: parzen x, q
  of HarrellDavis:
    if w.isNil:
      var w: seq[F]; harrellDavis w, x.len, q
      dot x[0].addr, w[0].addr, x.len
    else:
      dot x[0].addr, w, x.len
  of HodgesLehmann: hodgesLehmann x, q
  of HDW: hDW x, q

when isMainModule:
  when not declared(stdin): import std/[syncio, formatfloat]
  import std/[strutils, random], cligen
  when defined release: randomize()
  proc qtl(`method`=Parzen, cbα=0.05, ps: seq[float]) =
    ## Read column of numbers on stdin; Emit various quantiles for probabilities
    ## `ps`. Eg.: *echo 1 1 2 4|tr ' ' '\\n'|qtl 0 .5 1* writes *1.0 1.66 4.0*.
    ##
    ## If one `p`,integral>=2 (eg. ``210``), instead emit THAT many re-samples
    ## to stdout.  Eg. `|qtl 210|split -nl/21` creates 10 files each w/21 data.
    var x: seq[float]
    for line in stdin.lines: x.add parseFloat(line.strip)
    rsort x
    if ps.len == 1 and ps[0] == ps[0].int.float and ps[0] >= 2:
      for i in 1 .. int(ps[0]): echo x.quantile(rand(1.0))
    else:
      for i, p in ps:
        if i > 0: stdout.write " "
        stdout.write x.quantile(p, `method`, cbα)
      stdout.write '\n'
  include cligen/mergeCfgEnv; dispatch qtl, help={
    "method": """ p)arzen
 d(HarrellDavis)
 l(ParzenQmid-HodgesLehmann)
 b(both HD&HL, i.e. HD(Walsh)""",
    "cbα"   : "DKWM conf.band α"}
