import std/[strformat, strutils, sequtils, algorithm, math],
       spfun/[gauss, gamma], basicLA, svdx, min1d
template sum(F, i, it, expr): untyped = (var s: F; for i in it: s += expr; s)
type
  F8     = float
  CrVal* = enum xvGCV = "GCV", xvLOO = "LOO"

proc colDel(x: pointer; n, m, J, sz: int) = # In-place column deletion
  if m < 2: return                          # 2*n memmoves copying most data
  var X = cast[cstring](x)                  # One place that col-major sucks.
  let M = m - 1
  for i in 0 ..< n:
    if J > 0: moveMem X[M*i*sz].addr, X[m*i*sz].addr, J*sz
    if J < M: moveMem X[(M*i + J)*sz].addr, X[(m*i + J + 1)*sz].addr, (M - J)*sz

type CrossVal = object
  n: int        # number of rows of U/data points
  m: int        # number of predictors/params; dim of z,w
  z: seq[F8]    # z[j] sum(U_ij*y_i); Data vector in canonical coordinates
  u: seq[F8]    # u[i+n*j] ith element of jth left singular vector
  y: seq[F8]    # y[i] ith element of response vector
  s: seq[F8]    # s[j] singular values
  rss0: F8      # residual sum of squares for full model (input)
  df: F8        # effective degrees of freedom (output)

proc loo(o: var CrossVal; L: F8): F8 =
  # Return Leave-One-Out Cross-Validation score also known as "Allen's PRESS" as
  # a function of ridge parameter L & in-hand Singular Value Decomposition.
  let n = o.n; let m = o.m
  for j in 0 ..< m: o.z[j] = sum(F8, i, 0..<n, o.u[i + n*j]*o.y[i]) # zj = Uj*y
  var ss: F8
  for i in 0 ..< n:
    var Hii, yEst: F8
    for j in 0 ..< m:
      let adj = 1.0 / (1.0 + L / (o.s[j]*o.s[j]))
      yEst += o.u[i + n*j]*o.z[j]*adj
      Hii  += o.u[i + n*j]*o.u[i + n*j]*adj
    let pr = (o.y[i] - yEst)/(1.0 - Hii)            # prediction residual
    ss += pr*pr                                     # sum of its square
  return ss / n.F8                                  # LOO aka Allen PRESS

proc gcv(o: var CrossVal; L: F8): F8 =
  # Return Generalized Cross-Validation score as a function of ridge parameter L
  # &in-hand Singular Value Decomposition; Implementation follows Golub79 Eq2.3.
  let n = o.n; let m = o.m
  var df  = F8(n - m)
  var rss = 0.0
  if o.rss0 == 0.0:
    for j in 0 ..< m: o.z[j] = sum(F8, i, 0..<n, o.u[i + n*j]*o.y[i]) # zj=Uj*y
    for i in 0 ..< n:                               # r = y - yEst = y - U.z
      let yEst = sum(F8, j, 0..<m, o.u[i + n*j]*o.z[j])
      rss += (o.y[i] - yEst)*(o.y[i] - yEst)        # rss against OLS yEst
    o.rss0 = rss
  else: rss = o.rss0
  for j in 0 ..< m:
    let adj = 1.0 / (1.0 + o.s[j]*o.s[j] / L)       # L-adjustment for s[j]
    df  += adj
    rss += (adj*o.z[j])*(adj*o.z[j])
  o.df = df
  return n.F8*rss/(df*df)   # Adding this to rss & re-eval makes =~ LOOCV

proc linFit*[F](y, x: openArray[F]; n, m: int; b, u, s, v, r, h: var seq[F];
            thr: var F=1e-6; xv=xvLOO; log: File=nil): (F,F,F) {.discardable.} =
  ## Find best fit b in Xb=y; Caller pre-adjusts design matrix X&y for y weight.
  ## X is column major for dot product speed in SVD. Ie. x[i + n*j] indexes j-th
  ## predictor at i-th data point (&U similarly). Fills b[m] w/coefs, s with adj
  ## sing.value reciprocals. v.len>0 => fill w/Cov(b) estimate. r.len>0 => fill
  ## w/resids.  h.len>0 => fill w/n DIAG vals of Hat Matrix X*inv(X'*X)*X'. thr
  ## >0 is SV CLIP as *frac* of maxSV; <0 => ridge w/lambda -thr; ==0 => LOO|GCV
  ## -optimal ridge.  Returns (ssR, df, ssY) if r.len>0 else (0, df, 0).
  let doC = v.len > 0; let doR = r.len > 0; let doH = h.len > 0
  b.setLen m                                    # maybe alloc all needed RAM
  if u[0].addr != x[0].unsafeAddr:              # setLen & copy
    u.setLen x.len
    copyMem u[0].addr, x[0].unsafeAddr, x.len*F.sizeof
  s.setLen m
  v.setLen m*m                                  # need if need `b` => always
  if doR: r.setLen n
  if doH: h.setLen n
  let svMx = svdx(u, s, v, n, m)                # SVD LHS of X.b=y design eq
  if doH:
    for i in 0 ..< n: h[i] = sum(F, j, 0..<m, u[i + n*j]*u[i + n*j])
  var df = F(n - m)
  if thr > 0:                                   # Manual singular value clip
    for j in 0 ..< m: s[j] = if s[j] > thr*svMx: 1.0 / s[j] else: 0.0
  else:
    if thr < 0: thr = -thr                      # Manual ridge parameter
    else:                                       # GCV|LOO-optimal ridge param
      var o = CrossVal(n: n, m: m, z: newSeq[F](m), u: u, y: y.toSeq, s: s)
      proc xvScore(lam: F8): F8 = (if xv==xvGCV: gcv(o,lam) else: loo(o,lam))
      minBrent(thr, F8(0), F8(svMx)*F8(n), xvScore, tol=0.001, itMx=99)
    df = n.F
    for j in 0 ..< m:                           # adjust SVs; eff.deg.freedom
      var recipDen = 1.0 / (s[j]*s[j] + thr)
      df -= s[j]*s[j]*recipDen                  #NOTE: s = s / (s^2 + thr)
      s[j] = s[j]*recipDen    # Ridge adj reduces in lim to s=1/s w/thr>0 above
    if not log.isNil:
      log.write &"sThr: {thr:.4g} n-df: {n.F - df:.5g} 1/sAdj: "
      for j in 0 ..< m: log.write &"{s[j]:.4g}%s", (if j < m-1: " " else: "\n")
  for k in 0 ..< m:                             # Calc best fit params b[]
    b[k] = F(0)
    for i in 0 ..< n: b[k] += sum(F, j, 0..<m, v[k + m*j]*s[j]*u[i + n*j]*y[i])
  var ssR, ssY: F
  if doR or doC:                                # residuals if requested|needed
    let (_, vY) = mvar(y); ssY = vY*F(n)
    for i in 0 ..< n:                           # yEst: estimated y from x&b
      let yEst = sum(F, j, 0..<m, b[j]*x[i + n*j])
      if doR: r[i] = y[i] - yEst
      ssR += (y[i] - yEst)*(y[i] - yEst)
  if doC:                                       # Cov(b) = inv(X'X) = V.W^2.Vt
    var cov = newSeq[F](m*m)
    let redCsq = if df > 0.0: ssR/df else: 0.0  # reduced Chi-square
    for i in 0 ..< m:
      for j in 0 .. i:
        let Cij = redCsq * sum(F, k, 0..<m, v[i + m*k]*v[j + m*k] * s[k]*s[k])
        cov[m*i + j] = Cij; cov[m*j + i] = cov[m*i + j]
    v = cov
  (ssR, df, ssY)

proc normalize*[F](x: var openArray[F]; n,M: int; xfm: string; o,s: var seq[F])=
  ## Here `xfm[j]` specifies 'c': centering only; 'z': z-score (0mean,1var);
  ## 'm': min-max (data range->[0,1]).  For explanations see:
  ## www.listendata.com/2017/04/how-to-standardize-variable-in-regression.html.
  for j, ch in xfm:
    case ch.toLowerAscii
    of 'c':
      let (avg, _) = mvar(x[j*n ..< (j+1)*n])
      o[j] += avg
      vadd x[j*n].addr, n, -avg
    of 'z':
      let (avg, vr) = mvar(x[j*n ..< (j+1)*n])
      o[j] += avg; s[j] *= sqrt(vr)
      vadd x[j*n].addr, n, -avg
      vmul x[j*n].addr, n, F(1)/sqrt(vr)
    of 'm':
      let (mn, mx) = mnmx(x[j*n ..< (j+1)*n])
      o[j] += mn; s[j] *= mx - mn
      vadd x[j*n].addr, n, -mn
      vmul x[j*n].addr, n, F(1)/(mx - mn)
    else: discard

proc linFit*[F](X: var openArray[F]; n,M: int; b,u,s,v, r,h, oX,sX: var seq[F];
       xfm=""; thr: var F=1e-6; xv=xvLOO; log: File=nil):(F,F,F){.discardable.}=
  ## This wraps 2-input `linFit` to normalize predictors/response according to
  ## `xfm` returning the offset,scale in oX,sX.  Other parameters & return are
  ## the same, but this 1-input `X` variant assumes the first column of cap `X`
  ## is `y` (for convenience with non-flat pre-weighted input & naming offScales
  ## for `y`). `xfm` is described in `normalize` doc above.
  normalize(X, n, M, xfm, oX, sX)
  linFit(X[0..<n], X[n..^1], n,M-1, b,u,s,v, r,h, thr, xv, log)
         #^y slice, ^preds slice; Col-major makes this easy

proc linFit*[F](X: var openArray[F]; n:var int; M:int; b,u,s,v,r,h,oX,sX: var seq[F];
       trim=F(0), its=1; xfm=""; thr: var F=1e-6; xv=xvLOO;
       log: File=nil): (F,F,F) {.discardable.} =
  ## This wraps 1-input `linFit` to trim outliers adding `trim` threshold (in
  ## "N(0,1)-sdev 'units'") & `its` max iterations.  (Doing more than 2..4 its
  ## is probably ill-advised, but it is also not obvious 1 is always enough.)
  if trim > F(0): r.setLen n            # Csq result MUST be meaningful to trim
  for it in 0..its:
    result = linFit(X, n,M, b,u,s,v, r,h, oX,sX,xfm, thr, xv, log)
    if trim > F(0):
      let rtReducedCsq = sqrt(result[0]/result[1])
      var t = -gauss.qtl(trim/F(2)/F(n))*rtReducedCsq # E[outlier] => r units
      let n0 = n
      for i in 0 ..< n:
        if r[i].abs > t:
          colDel X[0].addr, M, n, i - (n0 - n), F.sizeof
          dec n
      if n == n0: break                 # No points trimmed => done
      if log != nil: log.write &"{n0-n} / {n0} > {t/rtReducedCsq} sds\n"

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
    if ixX[j]==0: result.add &"{b[j-1]:.9g}"; result.add sep # No *1. 4intercept
    else:
      result.add &"{b[j-1]:.9g} *"
      result.add fmtBasis(cols[j][0], ixX[j], o[j], s[j], sep)

proc Q*[F](df, ssR: F): F = F(1) - gammaI(df/F(2), ssR/F(2))
