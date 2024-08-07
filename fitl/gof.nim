## Implement metrics from the book Goodness-of-fit Techniques, 1986 edited by
## Ralph B. D'Agostino and Michael A. Stephens (referred to now as ds86).
import std/[math,random,algorithm,sugar],spfun/[gauss,gamma,binom],basicLA,ksamp
type              # Open Topology              Circular Topology
  GoFTest* = enum gfD  = "kolmogorovSmirnovD", gfV  = "vKuiper" , # L_infinity
                  gfW2 = "cramerVonMisesW2"  , gfU2 = "watsonU2", # L2 norm
                  gfA2 = "andersonDarlingA2"

  GoFMod* = enum gfFin = "finiteN", gfEst = "estimates"

const gofName*: array[GoFTest, string] = [ "D", "V", "W^2", "U^2", "A^2" ]
const noMod: set[GoFMod] = {}

func kolmogorovSmirnovPM[F](ps: seq[F]): (F, F) =
  let nInv = F(1.0/float(ps.len))                               # ds86 Eq 4.2
  for i, p in ps:       #NOTE: formula converted to 0-origin indexes.
    result[0] = max(result[0], F(i+1)*nInv - p)
    result[1] = max(result[1], p - F(i)*nInv)

func kolmogorovSmirnov*[F](ps: seq[F], mods=noMod): F =
  ## Kolmogorov-Smirnov max(D+,D-) modified by ds86 for estimated mean, var,
  ## finite n. `u01ize` first!
  let (dP, dM) = ps.kolmogorovSmirnovPM                         # ds86 Eq 4.2
  result = max(dP, dM); let n = ps.len.F
  if gfFin in mods: result *= sqrt(n) + F(0.12) + F(0.11)/sqrt(n)
  if gfEst in mods: result *= sqrt(n) - F(0.01) + F(0.85)/sqrt(n)

func kuiperV*[F](ps: seq[F], mods=noMod): F =
  ## Kuiper V=Dp+Dm modified by ds86 for est.mean,var,finite n. `u01ize` first!
  let (dP, dM) = ps.kolmogorovSmirnovPM
  result = dP + dM; let n = ps.len.F                           # ds86 Eq 4.2
  if gfFin in mods: result *= sqrt(n) + F(0.155) + F(0.24)/sqrt(n)
  if gfEst in mods: result *= sqrt(n) + F(0.050) + F(0.82)/sqrt(n)

func cramerVonMises*[F](ps: seq[F], mods=noMod): F =
  ## W^2 modified by ds86 for estimated mean,var, finite n. `u01ize` first!
  let n = ps.len
  let nI2 = F(0.5)/F(n) #NOTE: formula converted to 0-origin indexes.
  result = nI2/F(6) + sum0(j, n, (ps[j] - F(2*j + 1)*nI2)^2)    # ds86 Eq 4.2
  if gfFin in mods:
    result *= (F(1) - F(0.4)/F(n) + F(0.6)/F(n*n))*(F(1) + F(1)/F(n))
  # Stephens1970, ds86 Table 4.2 & Eq6.19 all say '*',but '/' in ^^^ E.g. 4.4.1
  if gfEst in mods: result *= F(1) + nI2

func watsonU2*[F](ps: seq[F], mods=noMod): F =
  ## Watson U^2 modified by ds86 for estim.mean,var, finite n. `u01ize` first!
  let n = ps.len; let nF = n.F
  let nI2 = F(0.5)/nF   #NOTE: Formula converted to 0-origin indexes.
  let mn = sum0(j, n, F(ps[j]))/nF                              # ds86 Eq 6.18
  result = nI2/F(6) - nF*(mn - 0.5)^2 + sum0(j, n, (ps[j] - F(2*j + 1)*nI2)^2)
  if gfFin in mods:
    result *= (F(1) - F(0.1)/F(n) + F(0.1)/F(n*n))*(F(1) + F(0.8)/F(n))
  # Stephens1970, ds86 Table 4.2 & Eq6.19 all say '*',but '/' in ^^^ E.g. 4.4.1
  if gfEst in mods: result *= F(1) + nI2

func andersonDarling*[F](ps: seq[F], mods=noMod): F =
  ## A^2 modified by ds86 for estimated mean,var, finite n. `u01ize` first!
  let n = ps.len                                                # ds86 Eq 4.2
  let nInv = F(1)/F(n)  #NOTE: Formula converted to 0-origin indexes.
  result = F(-n) - nInv*sum0(j, n, F(2*j + 1)*ln(ps[j]*(F(1) - ps[n-1-j])))
  if gfEst in mods: result *= F(1) + F(0.75)*nInv + F(2.25)*nInv*nInv

func zscore*[F](xs: var seq[F]; mnVr: (F, F)) =
  ## Unitize a data sample `xs[i]` to its z-scores `(xs[i] - mn)/sqrt(vr))`.
  let scl = F(1)/sqrt(mnVr[1])
  for x in mitems(xs): x -= mnVr[0]; x *= scl

func pitz*[F](zs: var seq[F]) =
  ## z-scores `N(0,1)->U(0,1)` via Probability Integral Transform (PIT).
  for z in mitems(zs): z = gauss.cdf(z)

func u01ize*[F](xs: var seq[F], mnVr: (F, F)) =
  ## Convert into Z-scores, sort, and PIT-transform to be U(0,1)
  xs.zscore mnVr; xs.sort; xs.pitz

func gofStat*[F](sample: var seq[F], mnVr: (F, F), gof=gfA2, mods=noMod,
                 u01d=false): F =
  ## Calculate goodness-of-fit stat for `sample` (clobbering `sample`).
  if not u01d: sample.u01ize(mnVr)
  case gof
  of gfD : sample.kolmogorovSmirnov(mods)
  of gfV : sample.kuiperV(mods)
  of gfW2: sample.cramerVonMises(mods)
  of gfU2: sample.watsonU2(mods)
  of gfA2: sample.andersonDarling(mods)

func gofStat*[F](sample: seq[F], mnVr: (F, F), gof=gfA2, mods=noMod,
                 u01d=false): F =
  ## Calculate goodness-of-fit stat for `sample` (preserving `sample`).
  var sample = sample
  sample.gofStat mnVr, gof, mods, u01d

template gofDistT(F, cdf, mV, gof, mods): untyped =
  let (mu, sig) = (mV[0], sqrt(mV[1]))
  var sample = newSeq[F](n)
  cdf.setLen 0
  for k in 1..m:
    for i in 0..<n: sample[i] = gauss.qtl(rand(1.0))*F(sig) + F(mu)
    cdf.add sample.gofStat(sample.mvars, gof, mods)
  cdf.sort
proc gofDist*(cdf: var seq[float32], n: int, mV=(0f32,1f32), gof=gfA2,
              mods=noMod, m=5000) =
  ## Fill `cdf` w/CDF of GoF stat `gof` for sample size `n`. Slow for large n*m.
  gofDistT(float32, cdf, mV, gof, mods)
proc gofDist*(cdf: var seq[float64], n: int, mV=(0f64,1f64), gof=gfA2,
              mods=noMod, m=5000) =
  ## Fill `cdf` w/CDF of GoF stat `gof` for sample size `n`. Slow for large n*m.
  gofDistT(float64, cdf, mV, gof, mods)

func prob*[F](cdf: seq[F], st: F): F =
  ## Return `P(x <= st)` {p-value}
  cdf.lowerBound(st)/cdf.len #XXX Interpolate via some local cubic polynom fit.
  #NOTE Also, stats were adjusted by Stephens to be ~ independent of n.

template gofTestTmpl(F, ps, mV, g, mods, m): untyped =
  result[0] = ps.gofStat(mV, g, mods, true)
  var cdf: seq[float32]
  cdf.gofDist(ps.len, mV, g, mods, m)
  result[1] = 1.0 - cdf.prob(result[0])
proc gofTest*[float32](ps: seq[float32], mV=(0f32,1f32), g=gfA2, mods={gfEst},
                       m=5000): (float32,float32)=
  ## Do whole test (stat+prob) calculation from u01ize'd probabilities `ps`.
  gofTestTmpl(float64, ps, mV, g, mods, m)
proc gofTest*[float64](ps: seq[float64], mV=(0f64,1f64), g=gfA2, mods={gfEst},
                       m=5000): (float64,float64)=
  ## Do whole test (stat+prob) calculation from u01ize'd probabilities `ps`.
  gofTestTmpl(float64, ps, mV, g, mods, m)

proc abs[F](x: seq[F]): seq[F] =
  result.setLen x.len; for i,e in x:result[i] = e.abs
proc ms[F](x: seq[F]): (F, F) =
  var sx, sx2: F
  for e in x: sx += e; sx2 += e*e
  result[0] = sx/x.len.F
  result[1] = 1.0/sqrt(sx2/x.len.F - result[0]^2)
proc u[F](x: var seq[F]) =
  let (m, s) = x.ms
  for i,e in x: (x[i] -= m; x[i] *= s)
proc dot[F](x, y: seq[F]): F = (for i,e in x: result += e*y[i])

proc whiteNoise*[F](x: openArray[F], K=5): F =
  ## pValue for W)hite N)oise test of arxiv.org/abs/2203.10405; Simultaneous
  ## test for non-zero ACF (non-I) *OR* linear heteroscedasticity (non-ID).
  var s = 0.0; let n = x.len.F
  for k in 1..K:
    var xk0 = x[0..^k]; var xk0a = xk0.abs; xk0.u; xk0a.u
    var xk1 = x[k..^1]; var xk1a = xk1.abs; xk1.u; xk1a.u
    s += (dot(xk0 , xk1)^2 + dot(xk0a, xk1a)^2 +
          dot(xk0a, xk1)^2 + dot(xk0 , xk1a)^2)/(n - k.F)
  χ²c(4.0*K.F, n*(n + 2.0)*s)           # χ² Survival Function/Complementary F()

proc identDist*[F](x: openArray[F], alpha=0.05, ident=2..3): F =
  ## Return first failed pValue for multiple hypothesis test of ident.b shuffles
  ## of 1/ident.a fractions being all from an identical sampling distribution.
  var x = collect(for e in x: e)        # Check Homogeneity/Identical Distro
  var ns: seq[int]
  for j in 0..<ident.a - 1: ns.add x.len div ident.a
  ns.add x.len - ns.sum
  for trial in 1..ident.b:
    x.shuffle; let (_, midR) = adkSim(x, ns, m=1000)
    let (_, pHi) = initBinomP(midR, 1000).est(1 - alpha)
    if pHi < alpha/ident.b.float:       # Bonferroni
      return pHi                        # No need to keep trying
  1.0

when isMainModule:
  import std/strformat, cligen, fitl/covar
  when not declared(File): import std/formatfloat

  proc `$`(xs: seq[float]): string =
    result.add $xs[0]
    for x in xs[1..^1]: result.add " "; result.add $x

  type Emit=enum eZ="z", ePITz="PITz", eStat="stat", eDist="dist", eProb="prob"

  proc gof(sample: seq[float], gofs: seq[GoFTest] = @[], adj: set[GoFMod]={},
        emit={eStat}, m=5000, knownM=0.0, knownV=0.0, pval=0.05, cf=3,
        ident=2..3, ran=true): int =
    ## Tests for a Gaussian shape with known|estimated parameters.  E.g.:
    ##  *gof -g,=,k,v,c,w,a -ep 15 16 17 18 19 19 20 20 21 22 22 23 23 23 24 27*
    ## If `prob` is in `emit` then exit status doubles as a boolean test for
    ## significant departure from a Gaussian at the *alpha* level `pval`.
    if ran: randomize()
    var sample = sample                 #NOTE: Sample variance to match ds86
    let mV   = if knownV == 0.0: sample.mvars else: (knownM, knownV)
    let adj  = if knownV == 0.0 and adj.len == 0: {gfEst} else: adj
    let gofs = if gofs.len == 0: @[ gfA2 ] else: gofs
    sample.zscore mV
    var nSignificantDepartures = 0
    if (let wn = sample.whiteNoise; wn < pval): # Independ & Linearly homosced.
      inc nSignificantDepartures;echo "significant LinearDepend/Heterosced: ",wn
    if (let pHi = sample.identDist(pval, ident); pHi < pval/ident.b.float):
      inc nSignificantDepartures; echo "significant Non-Identicality: ",pHi
    if eZ in emit: echo "zScores: ", sample
    sample.sort; sample.pitz
    if ePITz in emit: echo "PITz: ", sample
    var st: float; var cdf: seq[float]
    for g in gofs:
      if (emit*{eStat,eProb}).len != 0: st = sample.gofStat(mV, g, adj, true)
      if (emit*{eDist,eProb}).len != 0: cdf.gofDist(sample.len, mV, g, adj, m)
      let gName = if adj.len == 0: gofName[g] else: "m" & gofName[g]
      if eStat in emit: echo &"{gName}: {st:.4g}"
      if eProb in emit:
        let p = 1.0 - cdf.prob(st)
        if p < pval: inc nSignificantDepartures
        echo &"P({gName}>val|Gauss): {p:.4g}"
      if eDist in emit: echo &"cdf({gName}): ", cdf
    nSignificantDepartures

  dispatch gof, positional="sample", short={"knownM":'M',"knownV":'V',}, help={
   "sample": "x_1 x_2 .. x_n data sample",
   "gofs":
     "kolmogorovSmirnovD cramerVonMisesW2 `andersonDarlingA2` vKuiper watsonU2",
   "adj"   : "adjust GoF stat for: **estimates** finiteN",
   "emit"  : "emits: z, PITz, stat, dist, prob",
   "m"     : "number of n-samples to estimate CDF",
   "knownM": "known *mean* Gaussian; kV=0 => **estimate**",
   "knownV": "known *var*  Gaussian; 0 => **estimate**",
   "pval"  : "exit status = number of `prob` < `pval`",
   "cf"    : "test serial autoCorrFunc up to this lag",
   "ident" : "test identically distributed w/`b` trials of `a`-way splits",
   "ran"   : "randomize() for sampling"}
