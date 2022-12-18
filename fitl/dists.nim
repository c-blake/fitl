## Module defines a common framework for continuous probability distributions.
when not declared(stderr): import std/[syncio, formatfloat]
import std/[math, algorithm, random, critbits, strutils], spfun/[gauss, cauchy]
type
  T = float64 # Things may become generic over this; T to not confuse with F(x)
  CDist* = tuple[                 ## A 1-Dimensional Probability Distribution
    pdf, cdf, qtl: proc(x: T): T, ## Density, Cumulative Distro, qtl/inverse CDF
    gen    : proc(): T,           ## Sample pseudo-random deviates of this dist
    support: seq[T],              ## where non-zero (for plots, num.integ., etc.
    modes  : seq[T]]              ## Locations of modes (local maxima)

template newton(p0, cdf, pdf, support, x0, fTol) {.dirty.} = # ?df must use `x`!
  let p = p0          # Capture in case p0 is an expression like `rand(1.0)`
  var x = x0          # Newton's method should be an ok inverter/solver here.
  for it in 1..50:    # Quadratically cvgent. If 50 fails,likely nothing works
    let pr = (cdf) - p
    if abs(pr) < fTol and it > 1: return x
    x -= pr / (pdf)
  result = if p < T(0.5): support[0] else: support[^1]

proc mix*(compons: seq[(T, CDist, T, T)]; support: seq[T] = @[],
          modes: seq[T] = @[]): CDist =
  ## Make a distribution by mixing with component fractions (must total 1.0)
  ## various `CDist`s with location & scale shifts.  If empty, `support` &
  ## `modes` are inferred from components (with some assumptions, obviously).
  for (_, _, _, scale) in compons:
    if scale < 0: raise newException(ValueError, "negative scale")
  let pd = (proc(x: T): T =
    for (coef, dist, location, scale) in compons:
      result += coef / scale * dist.pdf((x - location)/scale))
  let cd = (proc(x: T): T =
    for (coef, dist, location, scale) in compons:
      result += coef * dist.cdf((x - location)/scale))
  var cumWt: seq[T] = @[compons[0][0]]  # This is captured in qtl closure
  for i in 1..compons.high: cumWt.add cumWt[^1] + compons[i][0]
  if abs(cumWt[^1] - T(1)) > 1e-6: stderr.write "unnormalized mixture!\n"
  result.gen = (proc(): T =
    let p = rand(1.0)
    for i, cw in cumWt:                 # Can binary search if MANY components.
      if p <= cw:                       # `i` is now the right component.
        let (coef, dist, location, scale) = compons[i]
        return location + scale*dist.qtl((cw - p)/coef)
    compons[^1][1].support[^1])         # p =~ 1.0 should be very rare
  let supp = if support.len > 0: support else: (
    var (supp0, supp1) = (T.high,T.low) # Low-side & High-side inferred support
    for (_, dist, location, scale) in compons:
      supp0 = min(supp0, dist.support[0]*scale  + location)
      supp1 = max(supp1, dist.support[^1]*scale + location)
    @[supp0, supp1])                    #NOTE: R code calls these "breaks"
  result.qtl = (proc(p: T): T = newton(p, cd(x), pd(x), supp, 0.11, 5e-6))
  result.pdf = pd; result.cdf = cd; result.support = supp
  result.modes = if modes.len > 0: modes else: (var modes: seq[T];
    for (coef, dist, location, scale) in compons:
      for mode in dist.modes:
        if (let m = mode*scale + location; m notin modes): modes.add m
    modes.sort; modes)

const z=T(0); const mH=T(-0.5); const pQ=T(0.25); const m6=T(-6);const p2=T(2)
const o=T(1); const pH=T(+0.5); const p4=T(+4)  ; const p6=T(+6);const oT=T(0.1)
const p20 = T(20); const o3=T(1)/T(3) #z)ero,o)ne,mX=-X, pX=+X; oT=(o)ne(T)enth

template PD(pd, cd, qt, supp, mo): untyped =    # For invertible formula distros
  ((proc(x {.inject.}: T): T {.closure.} = pd),
   (proc(x {.inject.}: T): T {.closure.} = cd),
   (proc(p {.inject.}: T): T {.closure.} = qt),
   (proc(): T {.closure.} = (let p {.inject.} = rand(1.0); qt)),
   supp, mo)

template PDn(pd, cd, supp, mo): untyped =       # For numerical invert distros
  ((proc(x {.inject.}: T): T {.closure.} = pd), # expr in terms of `x`
   (proc(x {.inject.}: T): T {.closure.} = cd), # expr in terms of `x`
   (proc(p: T): T {.closure.} = newton(p, cd, pd, supp, 0.125, 5e-6)),
   (proc(): T {.closure.} = newton(rand(1.0), cd, pd, supp, 0.125, 5e-6)),
   supp, mo)

# Define common `CDist`s;  Export so importers can do their own mixtures/etc.
let dU*: CDist = PD(if x < z: z elif x <= o: o else: z, # {}[0] needs a type
                    if x < z: z elif x <= o: x else: o,
                    if p < z: z elif p <= o: p else: z, @[z, o], @[pH])

let dN* = PD(gauss.pdf[T](x), gauss.cdf[T](x), gauss.qtl[T](p), @[m6, p6], @[z])

let dCauchy* = PD(cauchy.pdf[T](x), cauchy.cdf[T](x), cauchy.qtl[T](p),
                  @[-p20, p20], @[z])   # Instantiate generics from `spfun`

let dExp* = PD(if x < z: z else: exp(-x), if x < z: z else: o - exp(-x),
               if p < z: z else: -ln(o - p), @[z, T(10)], @[z])

let dLaplace* = PD(if x < z : pH*exp(x) else: pH*exp(-x),
                   if x < z : pH*exp(x) else: o - pH*exp(-x),
                   if p < pH: ln(p2*p)  else: -ln(p2 - p2*p), @[m6, p6], @[z])

let dTri* = PD(
  if x < -o: z elif x < z : o + x          elif x < o: o - x          else: z,
  if x < -o: z elif x < z : x + x*x*pH     elif x < o: x - x*x*pH     else: o,
  if p < z : z elif p < pH: o-sqrt(o-p2*p) elif p < o: sqrt(o+p2*p)-o else: o,
  @[-o, o], @[z])                       ## Triangular on [-1,1]

let dLogNormal* = PD(gauss.pdf[T](ln(x))/x, gauss.cdf[T](ln(x)),
                     exp(gauss.qtl[T](p)), @[z, p20], @[exp(-o)])

let distros* = {"U01":dU, "Exp": dExp,  #TODO Check supports&modes; Maybe plot?
  "Maxwell": PD(x*exp(mH*x*x), o-exp(mH*x*x), sqrt(-p2*ln(o-p)), @[m6,p6],@[o]),
  "Laplace": dLaplace,
  "Logistic": PD(exp(-x)/(o+exp(-x))^2,o/(1+exp(-x)),ln(p/(o-p)),@[m6,p6],@[z]),
  "Cauchy": dCauchy,
  "ExtVal"   : PD(exp(-exp(-x)-x), exp(-exp(-x)), -ln(-ln(p)), @[m6, p6], @[z]),
  "InfPeak"    : PD(o/sqrt(p4*x), x.sqrt, p*p, @[z, o], @[z]),
  "AsymPareto" : PD(pH/pow(x,T(1.5)), o-o/x.sqrt, o/(o-p)^2, @[o, p20], @[z]),
  "SymPareto"  : PD(pQ*pow(o+x.abs,T(-1.5)),
                    if x < 0: pH/sqrt(o+x.abs) else: o-pH/sqrt(o+x.abs),
                    if p<0.5: o-o/(o-abs(o-p2*p))^2
                    else    : o/(o-abs(p2*p-o))^2-o, @[-p20, p20], @[z]),
  "N01": dN, "logNormal": dLogNormal,
  "3BinHisto1" : mix(@[(pH,dU,mH,o), (pH,dU,T(-5),T(10))], modes = @[z]),
  "Matterhorn" : PD(o/(x.abs*ln(x.abs)^2),
    if x.abs < -exp(-p2): z elif x.abs<exp(-p2): pH - x.sgn.T/x.abs.ln else: o,
    (let P = p-pH; let PS = P.sgn.T; PS*exp(-PS/P)),
    @[-exp(-p2), -1e-3, +1e-3, exp(-p2)], @[z]),
  "Log"        : PDn(-ln(x), x - x*ln(x), @[z, o], @[z]), "Triangle": dTri,
  "Beta2,2"    : PDn(6*x*(1-x), (3-2*x)*x*x, @[z, o], @[pH]),
  "ChiSquare1" : PDn(exp(-pH*x)/sqrt(2*PI*x), erf(sqrt(pH*x)), @[z, p20], @[z]),
  "NormalCubed": PD(gauss.pdf[T](x.sgn.T*pow(x.abs, o3))*o3*pow(x^2,-o3),
                    gauss.cdf[T](x.sgn.T*pow(x.abs, o3)),
                    gauss.qtl[T](p)^3, @[mH, pH], @[z]),
  "InvExp"     : PD(exp(-o/sqrt(x))*pH*pow(x, T(-3)*pH), exp(-o/sqrt(x)),
                    o/ln(p)^2, @[z, p20], @[o/T(9)]),
  "Marronite"    : mix(@[(p2/T(3), dN,z,o), (o/T(3), dN,-p20,o/T(4))]),
  "SkewedBimodal": mix(@[(T(0.75), dN,z,o), (o/p4, dN,T(1.5),o/T(3))],
                     modes = @[0.0005447, 1.442525]), #e^-x^2/2+e^-9/2*(x-3/2)^2
  "Claw"         : mix(@[(pH, dN,z,o ), (oT, dN,-o,oT), (oT, dN,mH,oT),
                         (oT, dN,z,oT), (oT, dN,pH,oT), (oT, dN,o ,oT)]),
  "SmoothComb"   : mix(@[(T(32)/T(63), dN,T(-31)/T(21), T(32)/T(63)),
                         (T(16)/T(63), dN, T(17)/T(21), T(16)/T(63)),
                         ( T(8)/T(63), dN, T(41)/T(21),  T(8)/T(63)),
                         ( T(4)/T(63), dN, T(53)/T(21),  T(4)/T(63)),
                         ( T(2)/T(63), dN, T(59)/T(21),  T(2)/T(63)),
                         ( T(1)/T(63), dN, T(62)/T(21),  T(1)/T(63))]),
  "Caliper": PDn((let x2 = x.abs + (if x.abs < oT: oT else: z);
                 if x.abs>=oT and x.abs<=T(1.1): p2*(o-pow(x2-oT,o3)) else: z),
    if   x < -T(1.1): z
    elif x < -oT    : pH*(o-(p4*(x.abs-oT)-3*pow(abs(-x-oT), p4*o3)))
    elif x < oT     : pH
    elif x < T(1.1) : 0.5*(4*(x-0.1)-3*pow(abs(x-0.1), p4*o3))
    else: 1.0, @[-T(1.1), -oT, oT, T(1.1)], @[-oT, oT]),
  "TriModeU": mix(@[(pQ,dU,T(-20.1),oT), (pH,dU,-o,p2), (pQ,dU,p20,oT)]),
  "Sawtooth": mix(@[(oT, dTri,T(-9),o), (oT, dTri,T(-7),o), (oT, dTri,T(-5),o),
                    (oT, dTri,T(-3),o), (oT, dTri,T(-1),o), (oT, dTri,T(+1),o),
                    (oT, dTri,T(+3),o), (oT, dTri,T(+5),o), (oT, dTri,T(+7),o),
                    (oT, dTri,T(+9),o) ]),
  "BiLogPeak": PDn(mH*ln(abs(x*(o-x))),
                   pH*(-x*x.abs.ln + (o-x)*abs(o-x).ln) + x, @[z,o], @[z,o]),
  "Bimodal" : mix(@[(pH, dN,-o,o), (pH, dN,o,o)], modes = @[z]),
  "10Normal": mix(@[(oT, dN,T(-25),o), (oT, dN,T(-20),o), (oT, dN,T(-15),o),
                    (oT, dN,T(-10),o), (oT, dN,T( -5),o), (oT, dN,T( +0),o),
                    (oT, dN,T( +5),o), (oT, dN,T(+10),o), (oT, dN,T(+15),o),
                    (oT, dN,T(+20),o)]),
  "unif": PD(if x < -o: z elif x <= o: pH       else: z,
             if x < -o: z elif x <= o: pH*(x+o) else: o,
             if p < z: -o elif p <= o: p2*p-o   else: o, @[-o, o], @[z]),
  "epanechnikov": PDn(if x <= -o: z elif x<o: T(0.75)*(o - x*x)       else: z,
                      if x <= -o: z elif x<o: T(0.75)*x - T(0.25)*x^3 else: o,
                      @[-o, o], @[z]),
  "biweight": PDn(if x < -o: z elif x <= o: T(0.9375)*(o-x*x)^2 else: z,
    if x <= -o:z elif x<o:pH+x*(T(0.9375)+x*x*(T(0.1875)*x*x-T(0.625)))else:o,
                  @[-o, o], @[z]),
  "triweight": PDn(if x < -o: z elif x <= o: T(1.09375)*(o-x*x)^3 else: z,
    if x <= -o:z elif x<o:
      pH+x*((x*x*(T(0.65625)-T(0.15625)*x*x)-T(1.09375))*x*x+T(1.09375))else:o,
                   @[-o, o], @[z]),
  "cosine": PD(if x < -o: z elif x <= o: PI/p4*cos(PI*pH*x)     else: z,
               if x < -o: z elif x <= o: pH*(sin(PI*pH*x) + o)  else: o,
               if p < z: -o elif p <= o: - p2/PI*arcsin(o-p2*p) else: o,
               @[-o, o], @[z])
} ## A collection of distros for density estimation like R's `benchden` ordered
  ## like Berlinet,Devroye1994. Davies09 & Rozenholc09 both add more multimodal
  ## (some included here, others too poorly described).  Also, some KDE kernels.
# Davies09 has no non-plot defs of: Outlier, StronglySkewed, Trimodal, ExpMix,
# 8BinHisto, DiscreteComb, 24Normal. 8Bin looks ~irregular histo of Old Faithful
# data. Trimodal ~ Avg(Bimodal, N(0,1)) w/SOME ratio. DiscreteComb ~ SmoothComb,
# but 2*3 eq spaced modes w/avg seps.  24Normal..maybe 1/24 sep of 1/240sd?

proc match*(ds: openArray[(string, CDist)], pfx: string): int =
  ## Let end users use short names for any {} of CDist's, e.g. `distros`.
  try: return parseInt(pfx) - 1
  except CatchableError: discard
  let pfxNorm = pfx.toLowerAscii        # Normalized & Canonical strings
  var allNorm, allCanon, allNum: seq[string]
  var cbt: CritBitTree[(string, int)]
  for i, (name, _) in ds:
    allCanon.add name
    let norm = name.toLowerAscii
    if pfxNorm == norm: return i    # Use exact match even if prefix of another
    allNorm.add(norm)
    cbt[norm] = (name, i)
    allNum.add $(i+1) & " " & name
  var ks: seq[string]; var slot: seq[int]
  for v in cbt.valuesWithPrefix(pfxNorm):
    ks.add(v[0]); slot.add(v[1])
  if   ks.len == 1: return slot[0]
  elif pfx.len == 0:
    raise newException(IOError, "No CDist.  Choices:\n  " & allNum.join("\n  "))
  elif ks.len > 1: raise newException(IOError,
    ("Ambiguous prefix for continuous distro; \"$1\" matches:\n  $2\n" &
     "Empty string lists all") % [ pfx, ks.join("\n  ")])
  else: raise newException(IOError, ("No prefix match for \"$1\"" % [pfx]))

when isMainModule:      # A trivial command line driver mostly for testing.
  import cligen; when defined(release): randomize()
  proc dists(distro="U01", nSamp=0, test=false,plot=false,modes=false,v=false)=
    var dno = 0
    try: dno = distros.match(distro)
    except IOError as e: echo e.msg; quit(1)
    if v: stderr.write "Distro[", dno+1, "]: ", distros[dno][0], '\n'
    let dist = distros[dno][1]
    for i in 1..nSamp: echo dist.gen()
    if test:
      let x = 0.125     # All Ds != 0 here; A whole domain test is out of scope.
      let h = 5e-5      # for numerical derivative of CDF
      for (name, distro) in distros:
        let (pdf, cdf, qtl, _, _, _) = distro
        let p = pdf(x); let c = cdf(x); let q = qtl(c)
        let dc = (cdf(x + h) - c)/h
        if abs(p/dc - 1) > 3e-4:echo name," cdf|pdf mismatch; pd: ",p," dc: ",dc
        if abs(q  -  x) > 1e-4: echo name," cdf|qtl mismatch; cd: ",c," qt: ",q
    if plot:
      let scl = (dist.support[^1] - dist.support[0])/4096.0
      for i in 0..4095:
        let x = dist.support[0] + i.T * scl
        echo x, " ", if x > dist.support[0]: dist.pdf(x) else: 0.0
      echo dist.support[^1], " ", 0.0
    if modes: (for m in dist.modes: echo m)
  dispatch dists
