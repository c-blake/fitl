import std/[strformat, math, random, strutils, algorithm],
       cligen, cligen/[osUt, mslice, strUt],
       spfun/studentT, fitl/[basicLA, covar, linfit, gof, qtl]
when not declared(File): import std/[syncio, formatfloat]
type
  Gof* = enum gofR2="r2", gofXsq="xsq", gofPar="param", gofV="vKuiper",
              gofD="kolmogorovSmirnovD"               , gofU2="watsonU2",
              gofW2="cramerVonMisesW2", gofA2="andersonDarlingA2"
  Cov* = enum covLabel="label", covNorm="norm", covEst="est", covBoot="boot"
  F*   = float32

proc colCpy(dst, src: pointer; nD, nS, m, sz: int) =  # To make bootstrap data
  let dD  = cast[uint](nD * sz)                       #..in a col-major world.
  let dS  = cast[uint](nS * sz)
  var Dst = cast[uint](dst)
  var Src = cast[uint](src)
  let eod = Dst + dD * m.uint
  while Dst < eod:
    copyMem cast[pointer](Dst), cast[pointer](Src), sz
    Dst += dD; Src += dS

proc parseCols(cols: seq[string]): (string, seq[string]) =
  for col in cols:
    case col[0].toLowerAscii
    of 'z': result[0].add 'z'; result[1].add col[1..^1]
    of 'c': result[0].add 'c'; result[1].add col[1..^1]
    of 'm': result[0].add 'm'; result[1].add col[1..^1]
    else  : result[0].add '.'; result[1].add col

var iNm: string
var iFl: File
proc s(vs: openArray[MSlice], j: int): MSlice {.inline.} =
  vs[if j < 0: vs.len + j else: j - 1]

proc f(vs: openArray[MSlice]; i, j: int): F {.inline.} =
  try: result = if j == 0: F(1) else: F(parseFloat(s(vs, j)))
  except: stderr.write &"{iNm}:{i}:non-numeric column {j}: \"{$s(vs,j)}\"\n"

proc parseInp(cols: seq[string]; sep: Sep; X: var seq[F];
              ixW=0): (int, seq[int], string) =
  var xT: seq[F]        # This routine could be more efficient in SEVERAL ways.
  let (str, cols) = parseCols(cols); result[2] = str
  for col in cols: result[1].add parseInt(col.strip)
  let ixMx = max(ixW, result[1].max)
  let ixMn = min(ixW, result[1].min)
  let need = max(ixMx, -ixMn)
  var nums: seq[MSlice]
  var i = 0
  for (cs, n) in iFl.getDelims:
    inc i
    if n > 0 and cs[0] != '#':
      sep.split MSlice(mem: cs, len: n), nums
      if nums.len < need:
        stderr.write &"{iNm}:{i}:skipping too few columns ({nums.len}<{need})\n"
      else:
        let w = if ixW == 0: F(1) else: F(1)/nums.f(i, ixW)
        for j in result[1]: xT.add w*nums.f(i, j)
  result[0] = xT.len div cols.len
  X = xpose(xT, result[0], cols.len)

proc fmtCov*[F](s: string; v: seq[F]; m=0; norm=false, label=false): string =
  proc elt(i, j: int): F  =             # fmt cov/corr
    if norm:
      if i==j: sqrt(v[m*i + i])                       # std.errs
      else   : v[m*i + j]/sqrt(v[m*i + i]*v[m*j + j]) # corr.coefs
    else: v[m*i + j]
  result.add s & "-" & (if norm: "stderr-corr" else: "covariance") & " matrix\n"
  if label:                             # FULL SYMM. MATRIX WITH LABELS
    for i in 0..<m:
      for j in 0..<m: result.add &"{i} {j} {elt(i,j):#11.04g}\n"
  else:                                 # UNLABELED LOWER TRIANGULAR
    for i in 0..<m:
      result.add repeat(' ', 11*i); result.add &"{elt(i,i):#10.04g}"
      for j in i + 1..<m: result.add &" {elt(i,j):#+10.04g}"
      result.add '\n'

proc fmtBasis[F](ch: char; ix: int; o,s: F; sep: string): string =
  let ch = ch.toLowerAscii
  let pm = "+-"[int(o > F(0))]
  result.add (if   ch == 'c'      : &"(${ix} {pm} {abs(o)})"
              elif ch in {'z','m'}: &"(${ix} {pm} {abs(o)})/{s}"
              else: &"${ix}")
  result.add sep

proc fmtModel*[F](cols: seq[string]; ixX: seq[int]; M: int;
                  b, v, o, s: seq[F]): string =
  result.add fmtBasis(cols[0][0], ixX[0], o[0], s[0], "= ")
  for j in 1..<M:                       # Emit model: Y= & coeffs*Xj
    let sep = if j==M-1: "\n" else: " + "
    let bstr = fmtUncertainVal(b[j-1], sqrt(v[j-1 + (j-1)*(M-1)]), sigDigs=3)
    if ixX[j]==0: result.add bstr & sep # No *1. 4intercept
    else: result.add bstr & " *" & fmtBasis(cols[j][0], ixX[j], o[j], s[j], sep)

proc fmtPar[F](leading: string; bs, v, bT: seq[F]): string =
  result.add &"{leading}param(sigma)\tsig/par\t5% .. 95%\n"
  for j, b in bs:
    let sig = sqrt(v[bs.len*j + j])
    result.add &"{leading}{fmtUncertainMerged(b,sig)}\t{b/sig:.3g}\t"
    let digs = max(int(log(b/sig, 10.0) + 2.5), 2)
    if bT.len == 0:
      result.add formatFloat(b - 1.64485*sig, precision=digs); result.add " .. "
      result.add formatFloat(b + 1.64485*sig, precision=digs); result.add "\n"
    else:
      let n  = bT.len div bs.len
      var bb = bT[n*j ..< n*(j+1)]; bb.sort
      result.add formatFloat(bb.quantile(0.05),precision=digs);result.add " .. "
      result.add formatFloat(bb.quantile(0.95),precision=digs);result.add "\n"
  result.setLen result.len - 1

proc fitl*(cols: seq[string], file="-", delim="w", wtCol=0, sv=1e-8, xv=xvLOO,
           resids="", acf=0, boot=0, gof: set[Gof]={}, cov: set[Cov]={}, log="",
           trim=F(0), its=0) =
  ## Linear least squares parameter estimator for ASCII numbers in `file`.
  ## Default output is an awk/gnuplot-like formula w/best fit coefs to make ys
  ## from xs.  Options control extra output.  The input format is just:
  ##   <yCol> <basis1_column> [ <basis2_column> .. ]
  ##     y1      basis1(x1)   [    basis2(x1)   .. ] { Permutable cols.
  ##     y2      basis1(x2)   [    basis2(x2)   .. ]   1-origin numbers. }
  ##     .           .        [        .        .  ]
  ##     .           .        [        .        .  ]
  ## NOTE: Input ^$ blanks & lines beginning with '#' are skipped.  colNo '0' ->
  ## 1 for all data points -- useful as both flat reciprocal weight (aka sigma)
  ## & intercept col.  A colNo prefix of z => Z)SCORE it (mean0,var1); c => only
  ## C)ENTER it (mean0); m =>(x-min)/(max-min) {->[0,1]}.
  if cols.len < 2: raise newException(HelpError,"Too few columns; Full ${HELP}")
  let resF = if resids.len != 0: open(resids, fmWrite)  else: nil
  let logF = if log.len    != 0: open(log   , fmAppend) else: nil
  if file == "-": iNm = "stdin"; iFl = stdin
  else: iNm = file; iFl = open(file)
  let sep = initSep(delim)
  let M = cols.len                      # total num columns
  let m = M - 1                         # num x/predictor columns
  var X: seq[F]                                # Parse cols,text->Y,DesignMatrix
  var (n,ixX,xfm) = parseInp(cols,sep,X,wtCol) #..as well as centr&stdize ctrls
  if m > n: quit &"fewer data rows in {iNm} than columns", 3
  var o = newSeq[F](M)                  # Offset/Origin for each column (0.0)
  var s = newSeq[F](M, 1.0)             # Scale to divide by to normalize data
  var b = newSeq[F](m) # Do not clobber X w/u;  Since X is col-major aka bck2bck
  var u = X            #..cols, first col can be "y" & rest is still 1 mem block
  var w = newSeq[F](m)                  # Singular values/recips
  var v = newSeq[F](m*m)                # Right sing.vectors/Cov(b) matrix
  var r = newSeq[F](n)                  # Fit Residuals; Q: conditional alloc?
  var h: seq[F]                         #TODO Hat-matrix/influence functions
  var thr = F(sv)
  let (ssR,df,ssY) = linFit(X,n,M, b,u,w,v, r,h, o,s, trim,its,xfm, thr,xv,logF)
  echo fmtModel(cols, ixX, M, b, v, o, s)       # emit the model
  if resF != nil: (for i in 0..<n: resF.write &"{r[i]:.11g}\n")
  for lag in 1..acf:                            # Maybe emit resid AutoCorrFunc
    let ac = corrAuto(r, lag)
    echo &"ResidAutoCorr-Lag-{lag}: {ac:.5f}\t{F(1) - corrP(ac, n - lag):.5f}"
  if covEst in cov: echo fmtCov("estimated",v,m,covNorm in cov, covLabel in cov)
  if gofR2  in gof: echo &"r-squared: {F(1) - ssR/(ssY*s[0]):.6g}"
  if gofXsq in gof: echo &"Chi-sqr: {ssR:.4g} nu: {df:.4g} p: {Q(df, ssR):.5g}"
  var mV: (F, F)
  if ({gofD,gofW2,gofA2,gofV,gofU2}*gof).len>0: # *some* EDF-based residual test
    mV = r.mvars; r.u01ize mV                   # mean-vars => PITz just once
  proc fmtGf(nm: string; sP: (F, F); sDig,pDig: int): string =
    result.add nm;result.add ": "; result.add formatFloat(sP[0], precision=sDig)
    result.add "  "; result.add formatFloat(sP[1], precision=pDig)
    if sP[1] < 0.01: result.add '*'             # Q: take thresh as a param?
  if gofD  in gof: echo fmtGf("KSgaussRes" , r.gofTest(mV, gfD ), 4, 3)
  if gofW2 in gof: echo fmtGf("CvMgaussRes", r.gofTest(mV, gfW2), 4, 3)
  if gofA2 in gof: echo fmtGf("ADgaussRes" , r.gofTest(mV, gfA2), 4, 3)
  if gofV  in gof: echo fmtGf("KuiGaussRes", r.gofTest(mV, gfV ), 4, 3)
  if gofU2 in gof: echo fmtGf("WatGaussRes", r.gofTest(mV, gfU2), 4, 3)
  template slot: untyped = int(F(n)*rand(F(1))) # rand(n-1) unbiased but slow
  var bT: seq[F]                                # Collect `b` for raw `Cov(b)`
  if boot > 0:                                  # Bootstrapped cov(parameters)
    randomize()                                 # Fit synthetic new data sets..
    var Xp = newSeq[F](n*M)                     #..of the *same size*.
    var bK: seq[F]; var b = newSeq[F](m)
    r.setLen 0; h.setLen 0; let nn = xfm.needNormalize
    for k in 1..boot:                           # Gen data,reset,get&save coefs
      for i in 0..<n: colCpy Xp[i].addr, X[slot()].addr, n, n, M, F.sizeof
      v.zero; thr=sv; if nn: o.zero; s.set 1.0                # reset
      linFit(Xp,n,M, b,u,w,v, r,h, o,s,xfm, thr,xv,logF)      # get best fit b
      bK.add b                                                # save b
    bT.setLen bK.len; bT.xpose(bK, boot, m)
    covMat(v, bT, boot, m)                     # Replace `v` w/boostrapped cov
    if covBoot in cov:echo fmtCov("bootstrap",v,m,covNorm in cov,covLabel in cov)
  if gofPar in gof: echo &"Param Significance Breakdown:\n", fmtPar("  ",b,v,bT)
  (if resF != nil: resF.close); (if logF != nil: logF.close)
  if iFl != stdin: iFl.close # stdin must have seen EOF; So close not so wrong.

when isMainModule: dispatch fitl, help = {
  "cols"  : "1-origin-yCol xCol.. 0=>all 1s; *[cs]=>Centr/Std",
  "file"  : "input file; \"-\" => stdin",
  "delim" : "`initSep` input delimiter; w=repeated whitespc",
  "wtCol" : "1-origin sigma aka inverse weight column",
  "sv"    : "regularize: >0 SV clip <0 -manualRidge ==0 CV",
  "xv"    : "auto-ridge cross-validation score: GCV LOO",
  "resids": "log residuals to this pathname",
  "acf"   : "emit resid serial AutoCorrFunc up to this lag",
  "boot"  : "num bootstraps for Cov(b); 0=>estimated Cov(b)",
  "gof" : """emit goodness of fit diagnostics:
  r2: R^2; xsq: Chi-Square {aka SSR},df,pvalue
  param: parameter significance breakdown
  GoF tests residuals are Gaussian:
    kolmogorovSmirnovD cramerVonMisesW2
    andersonDarlingA2 vKuiper watsonU2""",
  "cov"   : "emit Cov(b) with flags: est, norm, label, boot",
  "log"   : "log to this path (trimming,model selection,..)",
  "trim": """trim pnts>="Nqtl(x/(2n)) sdevs" from regr surf
 [x=num pts expected if resids REALLY Normal]""",
  "its"   : "max trimming itrs; < 0 => until fixed point." }
