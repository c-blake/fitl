import std/[strformat, math, random, strutils], spfun/studentT, cligen,
       cligen/[osUt, mslice, strUt], basicLA, covar, linfit
type
  Gof* = enum gofR2="r2", gofCsq="csq", gofPar="param"
  Cov* = enum covLabel="label", covNorm="norm", covEst="est", covBoot="boot"
  F*   = float

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

proc s(vs: openArray[MSlice], j: int): MSlice {.inline.} =
  vs[if j < 0: vs.len + j else: j - 1]

proc f(vs: openArray[MSlice]; i, j: int): F {.inline.} =
  try: result = if j == 0: F(1) else: F(parseFloat(s(vs, j)))
  except: stderr.write &"stdin:{i}:non-numeric column {j}: \"{$s(vs,j)}\"\n"

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
  for (cs, n) in stdin.getDelims:
    inc i
    if n > 0 and cs[0] != '#':
      sep.split MSlice(mem: cs, len: n), nums
      if nums.len < need:
        stderr.write &"stdin:{i}:skipping too few columns ({nums.len}<{need})\n"
      else:
        let w = if ixW == 0: F(1) else: F(1)/nums.f(i, ixW)
        for j in result[1]: xT.add w*nums.f(i, j)
  result[0] = xT.len div cols.len
  X = xpose(xT, result[0], cols.len)

proc fmtPar[F](leading: string; bs, v: seq[F]): string =
  result.add &"{leading}par(sig)\tsig/par\n"
  for j, b in bs:
    let sig = sqrt(v[bs.len*j + j])
    result.add &"{leading}{fmtUncertainMerged(b,sig)}\t{b/sig:.3g}\n"
  result.setLen result.len - 1

proc fitl*(cols: seq[string], wtCol=0, delim="w", sv=1e-8, xv=xvLOO, resids="",
  acf=0, boot=0, gof:set[Gof]={}, cov:set[Cov]={}, log="", trim=F(0), its=0) =
  ## Linear least squares parameter estimator for ASCII numbers on stdin.
  ## Default output is an awk/gnuplot-like formula w/best fit coefs to make ys
  ## from xs.  Options control extra output.  The input format is just:
  ##   <yCol> <basis1_column> [ <basis2_column> .. ]
  ##     y1      basis1(x1)   [    basis2(x1)   .. ] { Permutable cols.
  ##     y2      basis1(x2)   [    basis2(x2)   .. ]   This is just for
  ##     .           .        [        .        .  ]   cols=1 2 3 4.. }
  ##     .           .        [        .        .  ]
  ## NOTE: Input ^$ blanks & lines beginning with '#' are skipped.  colNo '0' ->
  ## 1 for all data points -- useful as both flat reciprocal weight (aka sigma)
  ## & intercept col.  A colNo prefix of z => Z)SCORE it (mean0,var1); c => only
  ## C)ENTER it (mean0); m =>(x-min)/(max-min) {->[0,1]}.
  if cols.len < 2: raise newException(HelpError,"Too few columns; Full ${HELP}")
  let resF = if resids.len != 0: open(resids, fmWrite)  else: nil
  let logF = if log.len    != 0: open(log   , fmAppend) else: nil
  let sep = initSep(delim)
  let M = cols.len                      # total num columns
  let m = M - 1                         # num x/predictor columns
  var X: seq[F]                                # Parse cols,text->Y,DesignMatrix
  var (n,ixX,xfm) = parseInp(cols,sep,X,wtCol) #..as well as centr&stdize ctrls
  if m > n: quit "fewer data rows on stdin than columns", 3
  var o = newSeq[F](M)                  # Offset/Origin for each column (0.0)
  var s = newSeq[F](M, 1.0)             # Scale to divide by to normalize data
  var b = newSeq[F](m) # Do not clobber X w/u;  Since X is col-major aka bck2bck
  var u = X            #..cols, first col can be "y" & rest is still 1 mem block
  var w = newSeq[F](m)                  # Singular values/recips
  var v = newSeq[F](m*m)                # Right sing.vectors/Cov(b) matrix
  var r = newSeq[F](n)                  # Fit Residuals; Q: conditional alloc?
  var h: seq[F]                         #TODO Hat-matrix/influence functions
  var thr = sv
  let (ssR,df,ssY) = linFit(X,n,M, b,u,w,v, r,h, o,s, trim,its,xfm, thr,xv,logF)
  echo fmtModel(cols, ixX, M, b, v, o, s)       # emit the model
  if resF != nil: (for i in 0..<n: resF.write &"{r[i]:.11g}\n")
  for lag in 1..acf:                            # Maybe emit resid AutoCorrFunc
    let ac = corrAuto(r, lag)
    echo &"ResidAutoCorr-Lag-{lag}: {ac:.5f}\t{F(1) - corrP(ac, n - lag):.5f}"
  if covEst in cov: echo fmtCov("estimated",v,m,covNorm in cov, covLabel in cov)
  if gofR2  in gof: echo &"r-squared: {F(1) - ssR/(ssY*s[0]):.6g}"
  if gofCsq in gof: echo &"Chi-sqr: {ssR:.6g} nu: {df:.4g} p: {Q(df,ssR):.5g}"
  template slot: untyped = int(F(n)*rand(F(1.0))) # rand(n-1) unbiased but slow
  if boot > 0:                                  # Bootstrapped cov(parameters)
    randomize()                                 # Fit synthetic new data sets..
    var Xp = newSeq[F](n*M)                     #..of the *same size*.
    var bK = newSeq[F]()                        # Collect `b` for raw `Cov(b)`
    for k in 1..boot:
      for i in 0..<n: colCpy Xp[i].addr, X[slot()].addr, n, n, M, F.sizeof #Gen
      thr=sv; r.setLen 0; v.setLen 0; h.setLen 0; var u = newSeq[F](n*m) #Reset
      for j in 0..<M: o[j] = 0.0; s[j] = 1.0                  #offset&scale,too
      linFit(Xp,n,M, b,u,w,v, r,h, o,s, xfm, thr,xv,logF)               #Fit
      bK.add b                                                          #Record
    covMat(v, xpose(bK, boot, m), boot, m)      # Replace `v` w/boostrapped cov
    if covBoot in cov:echo fmtCov("bootstrap",v,m,covNorm in cov,covLabel in cov)
  if gofPar in gof: echo &"Param Significance Breakdown:\n", fmtPar("  ", b, v)
  if resF != nil: resF.close
  if logF != nil: logF.close

when isMainModule: dispatch fitl, help = {
  "cols"  : "1-origin-yCol xCol.. 0=>all 1s; *[cs]=>Centr/Std",
  "wtCol" : "1-origin sigma aka inverse weight column",
  "delim" : "`initSep` input delimiter; w=repeated whitespc",
  "sv"    : "regularize: >0 SV clip <0 -manualRidge ==0 CV",
  "xv"    : "auto-ridge cross-validation score: GCV LOO",
  "resids": "log residuals to this pathname",
  "acf"   : "emit resid serial AutoCorrFunc up to this lag",
  "boot"  : "num bootstraps for Cov(b); 0=>estimated Cov(b)",
  "gof" : """emit goodness of fit diagnostics:
  r2: R^2; csq: Chi-Square {aka SSR},df,pvalue
  param: parameter significance breakdown""",
  "cov"   : "emit Cov(b) with flags: est, norm, label",
  "log"   : "log to this path (trimming,model selection,..)",
  "trim": """trim pnts>="Nqtl(x/(2n)) sdevs" from regr surf
 [x=num pts expected if resids REALLY Normal]""",
  "its"   : "max trimming itrs; < 0 => until fixed point." }
