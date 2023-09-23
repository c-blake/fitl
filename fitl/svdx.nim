## Not so rotten self-contained SVD in 42 lines of non-comment/blank code
from math import sqrt, copySign; from basicLA import dot, sum0

proc jacobi[F](a, b, c: F): (F, F) {.inline.} =
  if abs(2.0 * c) != 0.0:               # return (cos,sin)(Jacobi Angle)
    let z = (b - a) / (2.0 * c)
    let t = copySign(1.0, z) / (abs(z) + sqrt(1.0 + z * z))
    result[0] = 1.0 / sqrt(1.0 + t * t)
    result[1] = -result[0] * t
  else:
    result[0] = 1.0 # result[1] = 0.0   # cos=1, sin=0; cos^2+sin^2=1

proc rot2[F](x, y: var F; co, sn: F) {.inline.} =
  let (u, v) = (x, y)                   # Apply 2D Givens Rotation
  x =  co*u + sn*v
  y = -sn*u + co*v

template xe(i, j: untyped): untyped = xu[i + n*j] # Elt accessor templates
template ve(i, j: untyped): untyped = v[i  + m*j]

proc svdx*[F](xu, s, v: var openArray[F]; n, m: int; tol=1e-6; mxIt=40): F =
  ## Thin/econ SVD (factor to u*s*vT) of a general N*M column major matrix xu.
  ## xu is replaced by U on return; Cols of U&V are Left&Right singular vectors.
  ## s & v must be presized to m and m*m, respectively.  v.len==0 skips calc of
  ## right singular vectors (e.g. for symm X).  tol=max off-diag matrix element
  ## (in units of sqrt(Si*Sj)) to converge.  Returns largest singular value.
  ## NOTE 1: Unsorted - sort[(s,j)..] to find the permutation if you need that.
  ## NOTE 2: u*(s>=0)*vT has 2^m sign choices; A flip in v cancels one in u.
  let doV = v.len > 0
  let t = tol * tol
  if doV:
    for i in 0 ..< m: ve(i, i) = 1.0    # V = I
  for it in 1..mxIt:
    var eMx = F(0)                      # Hestenes Algo which is just cyclic
    for i in 0 ..< m:                   #..Givens Rotations by Jacobi Angles.
      for j in i + 1 ..< m:             # std upper triangle
        let a = dot(xe(0, i).addr, xe(0, i).addr, n) # compute (Xt*xe)_ij
        let b = dot(xe(0, j).addr, xe(0, j).addr, n)
        let c = dot(xe(0, i).addr, xe(0, j).addr, n)
        let e = c * c / (a * b)         # a, b, c all >= 0
        eMx = max(eMx, e)
        if e > t:                       # Find angle to annihilate off diag elts
          let (co,sn) = jacobi(a, b, c) #..and then apply it all the way down u.
          for k in 0 ..< n  : rot2(xe(k, i), xe(k, j), co, sn)
          if doV:                       # optionally also apply to v cols
            for k in 0 ..< m: rot2(ve(k, i), ve(k, j), co, sn)
    if eMx < t: break                   # converged
  var sMax: F
  for j in 0 ..< m:                     # compute singular values
    s[j] = sqrt(dot(xe(0, j).addr, xe(0, j).addr, n))
    sMax = max(sMax, s[j])
    let sInv = 1.0 / s[j]               # also normalize U
    for i in 0 ..< n: xe(i, j) *= sInv
  sMax

template toOA(p, n) = toOpenArray[F](p, 0, n - 1)
proc polyFit*[F](b: var seq[F]; x,y,w: openArray[F], m=3, svThr=1e-6, tol=1e-6)=
  ## Fill `b` with `m` least-squares-best fit coefs in ascending order by power
  ## of polynomial model for y_i(x_i) maybe weighted by w_i.  w.len==0 => w_i=1.
  ## `svThr` is a singular value threshold in units of the max SV.  Smaller SVs
  ## are set to 0 for numerical stability even with (near) collinearities.
  let n = x.len
  if n < m: raise newException(ValueError, "n < m")
  if y.len != n: raise newException(ValueError, "y.len != x.len")
  if w.len>0 and w.len!=n: raise newException(ValueError, "w.len != x.len")
  if m < 1: raise newException(ValueError, "num.coefs m < 1")
  b.setLen (1 + 1 + m)*m + n*(m + 1)    # Output b, s, v, weighted Design Matrix
  var xuy = cast[ptr UncheckedArray[F]](b[(m + 2)*m].addr)
  for i in 0..<n:                       # Populate Design Matrix
    var val = if w.len > 0: w[i] else: F(1)
    xuy[i + n*0] = val                  # 0-th power of x[i] is 1.0
    for j in 1..<m:
      val *= x[i]; xuy[i + n*j] = val   # Build up & assign m-th power of x[i]
    xuy[i + n*m] = y[i]*(if w.len > 0: w[i] else: F(1)) # weight y[i]
  var s = cast[ptr UncheckedArray[F]](b[m].addr)   # skip b
  var v = cast[ptr UncheckedArray[F]](b[m+m].addr) # skip b & s
  let thr = svThr*svdx(toOA(xuy, n*m), toOA(s, m), toOA(v, m*m), n, m, tol)
  for j in 0..<m: s[j] = if s[j] < thr: F(0) else: F(1)/s[j]
  for k in 0..<m:                       # Calc best fit coefs (often beta)
    b[k] = sum0(j,m, sum0(i,n, v[k + m*j]*s[j]*xuy[i + n*j]*xuy[i + n*m]))
  b.setLen m

when isMainModule:
  when not declared(stdout): import std/[syncio, formatfloat]
  proc fmt[F](x: openArray[F]; n, m: int): string =
    for i in 0 ..< n:                   # format column major matrix
      for j in 0 ..< m: result.add ' '; result.add $x[i + n*j]
      result.add '\n'
  let x = [ 1.0f32,2,3,4, 5,6,7,8, 9,10,11,12 ] # 4rows(3 cols) COLUMN MAJOR
  var u = x                             # Make a copy so we can check later
  var s = newSeq[float32](3)
  var v = newSeq[float32](9)
  discard svdx(u, s, v, 4, 3)
  stdout.write "u cols Lv:\n",u.fmt(4,3),"s: ",s,"\n","v cols Rv:\n",v.fmt(3,3)
  from basicLA import sum0
  var e = 0f32      # U.S.Vt = sum_kl U_ik*S_kl*Vt_lj = U_ik*D_kl*s_k*Vt_jl
  for i in 0..<4:   #        = sum_k s_k*U_ik*V_jk
    for j in 0..<3: e=max(e, abs(x[i+4*j] - sum0(k,3, s[k]*u[i+4*k]*v[j+3*k])))
  echo "max|X-U.S.Vt|: ", e
