## Not so rotten self-contained SVD in 42 lines of non-comment/blank code
from math import sqrt, copySign; from basicLA import dot

proc jacobi[F](a, b, c: F): (F, F) {.inline.} =
  if abs(2.0 * c) != 0.0:               # return (cos,sin)(Jacobi Angle)
    let z = (b - a) / (2.0 * c)
    let t = copySign(1.0, z) / (abs(z) + sqrt(1.0 + z * z))
    result[0] = 1.0 / sqrt(1.0 + t * t)
    result[1] = -result[0] * t
  else:
    result[0] = 1.0; result[0] = 1.0

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
  let doV = v.len > 0
  let t = tol * tol
  if doV:
    for i in 0 ..< m: ve(i, i) = 1.0    # V = I
  for it in 1..mxIt:
    var eMx = F(0)                      # Hestenes algo which is just cyclic
    for i in 0 ..< m:                   #..Givens Rotations by Jacobi Angles.
      for j in i + 1 ..< m:             # std upper triangle
        let a = dot(xe(0, i).addr, xe(0, i).addr, n) # compute (Xt*xe)_ij
        let b = dot(xe(0, j).addr, xe(0, j).addr, n)
        let c = dot(xe(0, i).addr, xe(0, j).addr, n)
        let e = c * c / abs(a * b)
        eMx = max(eMx, e)
        if e > t:                       # calc Jacobi rot; apply to X
          let (co, sn) = jacobi(a, b, c)
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

when isMainModule:
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
