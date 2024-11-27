import std/[math, random, algorithm, stats], spfun/[gauss, digamma]
type F = float32            # Can make [F] a common generic parameter someday
# Common params: `x`: input; `n`: #pairs; `nB`: #boxes; `scl`: scale factor;
# `eGinv`: box size; `c(1|2)`: 1st|2nd component for making 2-D grid; `box`:
# each array elt is num of last point in box; `lis`: each array elt is num of
# prev points in box|-1; `mxi`: cumulative n(pnts) in box; eFoo =~ epsilon_Foo

func ssDup[T](n, d2, o: int, src: seq[seq[T]]): seq[seq[T]] = # seq[seq[]] dup
  result.setLen d2
  for d in 0 ..< d2:
    result[d] = newSeq[T](n)
    copyMem result[d][0].addr, src[o + d][0].addr, n*T.sizeof

proc xcol(x: seq[seq[F]]; i, dim: int): seq[F] =
  result.setLen dim; for d in 0 ..< dim: result[d] = x[d][i]

func make_box1(x: seq[F]; n: int; scl: F; nB: int; box,lis,mxi: var seq[int]) =
  for i in 0..nB: box[i] = -1; mxi[i] = 0   # Make 1-D box
  for i in 0 ..< n:
    let ix = int(x[i]*scl)
    lis[i] = box[ix]; box[ix] = i
    mxi[ix].inc
  for i in 1..nB: mxi[i] += mxi[i-1]

func make_box2ind(x: var seq[seq[F]]; dim, n, c1, c2, nB, eGinv: int;
                  ind, lis: var seq[int]; box: var seq[seq[int]]) =
  let ib = nB - 1 # Make 2-D box; Re-orders! but saves idx of orig data->`ind`.
  var xx = ssDup(n, dim, 0, x)
  for ix in 0 ..< nB: (for iy in 0..<nB: box[ix][iy] = -1)
  for i in 0 ..< n:
    let (ix,iy) = (int(x[c1][i]*eGinv.F) and ib, int(x[c2][i]*eGinv.F) and ib)
    lis[i] = box[ix][iy]; box[ix][iy] = i
  var i = -1
  for ix in 0 ..< nB:
    for iy in 0 ..< nB:
      var ixy = box[ix][iy]
      while ixy >= 0:
        inc i
        for d in 0 ..< dim: x[d][i] = xx[d][ixy]
        ind[ixy] = i; ixy = lis[ixy]
      box[ix][iy] = -1
  for i in 0 ..< n:
    let (ix,iy) = (int(x[c1][i]*eGinv.F) and ib, int(x[c2][i]*eGinv.F) and ib)
    lis[i] = box[ix][iy]; box[ix][iy] = i

func neiE1(x: seq[F]; scl: F; i,nB: int; eps: F; box,lis,mxi: seq[int]): int =
  let xc = x[i]                 # Count nbors of point i in eps-nborhood in 1-D
  let mp = min(nB, int((xc + eps)*scl))
  let mm = max(0, int((xc - eps)*scl))
  var mi = box[mp]
  while mi >= 0:
    let dd = x[mi] - xc
    if dd.abs <= eps: inc result
    mi = lis[mi]
  if mm >= mp: return result - 1
  mi = box[mm]
  while mi >= 0:
    let dd = xc - x[mi]
    if dd.abs <= eps: inc result
    mi = lis[mi]
  inc result, mxi[mp - 1] - mxi[mm]
  result - 1

iterator els(box: seq[seq[int]]; lis: seq[int]; ix2, iy1, ib: int): int =
  var el = box[ix2][iy1 and ib]
  while el != -1:
    yield el; el = lis[el]

template gridScan(eps, jj, ix, iy, ib, maxDscan) =
  while eps > eGrid*(jj.F - 1): # Outer loop logic for maxDscans in `nei[KE]`
    let step = if jj != 0: 2*jj else: 1
    for ix1 in ix - jj .. ix + jj:
      let ix2 = ix1 and ib
      for iy1 in countup(iy - jj, iy + jj, step): maxDscan ix2, iy1
    for ix1 in countup(ix - jj, ix + jj, step):
      let ix2 = ix1 and ib
      for iy1 in iy - jj + 1 .. iy + jj - 1: maxDscan ix2, iy1
    inc jj
    if jj == nB div 2: break
  if jj == nB div 2:            # Half of the layer
    for ix1 in ix - jj ..< ix + jj: maxDscan ix1 and ib, iy - jj
    let ix2 = (ix - jj) and ib
    for iy1 in iy - jj + 1 .. iy + jj - 1: maxDscan ix2, iy1

func neiK(x: seq[seq[F]]; dim, c1, c2, nB, i: int; eGrid: F; k: int;
          box: seq[seq[int]]; lis: seq[int]): seq[int] =
  let ib = nB - 1               # Search for k nbors of point i in dim-D ..
  let xx = xcol(x, i, dim)      # eGrid: size of grid; returns ixes(k nbors).
  let (ix, iy) = (int(xx[c1]/eGrid) and ib, int(xx[c2]/eGrid) and ib)
  var jj: int                   # Zero by default
  result.setLen k + 1
  var dn = newSeq[F](k + 1)
  for kk in 1..k: dn[kk] = F.high
  template maxDscan(ix2, iy1) =
    for el in els(box, lis, ix2, iy1, ib):
      if el != i:
        var dd = abs(xx[0] - x[0][el])
        for d in 1 ..< dim:
          if (let dy = abs(xx[d] - x[d][el]); dy > dd): dd = dy
        if dd < dn[k]:
          var kk = k
          while dd < dn[kk]:
            if kk < k: dn[kk + 1] = dn[kk]; result[kk + 1] = result[kk]
            dec kk
          dn[kk + 1] = dd; result[kk + 1] = el
  gridScan dn[k], jj, ix, iy, ib, maxDscan

func neiE(x: seq[seq[F]]; dim, c1, c2, nB, i: int; eGrid, eps: F;
          box: seq[seq[int]]; lis: seq[int]): int =
  let ib = nB - 1               # Count nbors in eps-nborhood of point i in
  let xx = xcol(x, i, dim)      # dim-D. eGrid: size of grid
  let (ix, iy) = (int(xx[c1]/eGrid) and ib, int(xx[c2]/eGrid) and ib)
  var jj, nx: int
  template maxDscan(ix2, iy1) =
    for el in els(box, lis, ix2, iy1, ib):
      var dd = abs(xx[0] - x[0][el])
      for d in 1 ..< dim:
        if (let dy = abs(xx[d] - x[d][el]); dy > dd):
          dd = dy
          if dd > eps: break
      if dd <= eps: inc nx
  gridScan eps, jj, ix, iy, ib, maxDscan
  nx - 1

proc mir_xnyn(x: var seq[seq[F]]; dimx, dimy, n, k: int; scl: seq[F]): F =
  let d2   = dimx + dimy        # Compute mutual info among vectors by rectangle
  let box1 = n - 5              # method for one k; IN: x: (dimx+dimy)-n-vectors
  var xc = newSeq[F](d2)        # dim(x|y): Dim of (x|y)vector; k: max num nbors
  var nB = 1; while nB*nB*k < 2*n: nB *= 2
  let eGinv = nB div 4; let eGrid = 1.0/eGinv.F
  var xx, yy: seq[seq[F]]
  var boxx, boxy, box: seq[seq[int]]
  var lisx,boxx1,lisx1,mxi, lisy,boxy1,lisy1,myi, lis,ind,indx,indy: seq[int]
  var scalx, scaly: F
  if dimx>1:
    xx = ssDup(n, dimx, 0, x)   # Save x to xx if data would be re-ordered
    boxx.setLen nB; for i in 0 ..< nB: boxx[i] = newSeq[int](nB)
    lisx.setLen n
  else: boxx1.setLen box1 + 1; lisx1.setLen n; mxi.setLen box1 + 1
  if dimy>1:
    yy = ssDup(n, dimy,dimx, x) # Save x to yy if data would be re-ordered
    boxy.setLen nB; for i in 0 ..< nB: boxy[i] = newSeq[int](nB)
    lisy.setLen n
  else: boxy1.setLen box1 + 1; lisy1.setLen n; myi.setLen box1 + 1
  box.setLen nB
  for i in 0 ..< nB: box[i] = newSeq[int](nB)
  lis.setLen n; ind.setLen n; indx.setLen n; indy.setLen n
  make_box2ind x, d2, n, 0, dimx, nB, eGinv, ind, lis, box
  if dimx>1: make_box2ind xx, dimx, n, 0, dimx-1, nB, eGinv, indx, lisx, boxx
  else: scalx = scl[0]; make_box1 x[0], n, scalx, box1, boxx1, lisx1, mxi
  if dimy>1: make_box2ind yy, dimy, n, 0, dimy-1, nB, eGinv, indy, lisy, boxy
  else: scaly = scl[dimx]; make_box1 x[dimx], n, scaly, box1, boxy1, lisy1, myi
  var dxy = 0.0
  for i in 0 ..< n:
    for d in 0 ..< d2: xc[d] = x[d][ind[i]]
    let nn = neiK(x, d2, 0, dimx, nB, ind[i], eGrid, k, box, lis)
    var ex = 0.0
    for d in 0 ..< dimx:
      for kk in 1..k: (if (let dx=abs(xc[d] - x[d][nn[kk]]); dx>ex): ex = dx)
    var ey = 0.0
    for d in dimx ..< d2:
      for kk in 1..k: (if (let dy=abs(xc[d] - x[d][nn[kk]]); dy>ey): ey = dy)
    let nx=if dimx==1: neiE1 x[0], scalx, ind[i], box1, ex, boxx1, lisx1, mxi
           else: neiE xx, dimx, 0, dimx-1, nB, indx[i], eGrid, ex, boxx, lisx
    let ny=if dimy==1: neiE1 x[dimx], scaly, ind[i], box1, ey, boxy1, lisy1, myi
           else: neiE yy, dimy, 0, dimy-1, nB, indy[i], eGrid, ey, boxy, lisy
    dxy += psi(nx) + psi(ny)    #echo "nx: ", nx, " ny: ", ny
  psi(n) + (psi(k) - 1.0/k.F) - dxy/n.F     # middle term = "phi[k]"

func zscoreScale(x: var seq[F]): F =
  var me, s: F; let n = x.len   # Make mean0,var1 & non-negative
  for i in 0 ..< n: me += x[i]  # This is what KSG 2004 recommends
  me /= n.F
  for i in 0 ..< n: s += (x[i] - me)^2
  s /= n.F - 1; s = 1.0/(sqrt s)
  var (mn, mx) = (F.high, F.low)
  for i in 0 ..< n:
    x[i] = s*(x[i] - me)
    if x[i] < mn: mn = x[i]
    if x[i] > mx: mx = x[i]
  for i in 0 ..< n: x[i] -= mn
  F(n - 5)/(mx - mn)            #echo "scl=", result

proc normalScoreScale*(pntsD: var seq[F]): F =
  let n  = pntsD.len            # This is what Holmes 2019 recommends.  Their..
  var xs = newSeq[(F, int)](n)  #.. "reparameterize" is aka "normal scoring".
  for i, x in pntsD: xs[i] = (x, i)
  xs.sort
  let q0 = gauss.qtl(0.5/n.F)   # Similar mean0,var1 but 0-offset as zscore
  for i in 0 ..< n: pntsD[xs[i][1]] = q0 + gauss.qtl((i.F + 0.5)/n.F)
  F(n - 5)/(gauss.qtl((n.F - 0.5)/n.F) - q0)

proc miKNN*(pnts: seq[F]; dimx, dimy, k: int; score=false): F =
  ## Mutual Info Estimate (bits) from data set of dimx-dimy vector pairs by k-th
  ## near-nbor rectangle method.  `pnts` is `n` back-to-back dimx+dimy vectors.
  let d2 = dimx + dimy
  let n = pnts.len div d2 # echo "n: ",n," k: ",k," dimx: ",dimx," dimy: ",dimy
  if pnts.len mod d2 != 0: raise newException(IOError,"data-dimx+dimy mismatch")
  var x = newSeq[seq[F]](d2)
  var scl = newSeq[F](d2)
  for d in 0 ..< d2:            # Set up whole coordinate columns
    x[d] = newSeq[F](n)
    for i in 0 ..< n: x[d][i] = pnts[d2*i + d]
    scl[d] = if score: x[d].normalScoreScale else: x[d].zscoreScale
  mir_xnyn(x, dimx, dimy, n, k, scl)/ln(2.0)    # Info is in bits not "nats"

iterator subsampled*(pnts: seq[F]; dimx, dimy, k, nS: int; score=false): F =
  ## Yield subsampled `miKNN` estimates for random 1/nS subsets of input `pnts`,
  ## sharing all other parameters.
  let n = pnts.len div nS   # Q: Insist `nS > 1`?
  var all = newSeq[int](pnts.len)
  for i in 0 ..< pnts.len: all[i] = i
  all.shuffle
  var ps = newSeq[F](n*(dimx + dimy))
  let bytes = (dimx + dimy)*F.sizeof
  for s in 1..nS:
    for i in 0 ..< n: copyMem ps[i].addr, pnts[all[(s - 1)*n + i]].addr, bytes
    yield miKNN(ps, dimx, dimy, k, score)

proc miErr*(pnts: seq[F]; dimx, dimy, k: int; score=false, err=4): (F, F) =
  result[0] = miKNN(pnts, dimx, dimy, k, score)
  var varSum = 0.0; var niSum = 0   # Paper derives assuming variance ~ χ².
  for ni in 2..err:                 # Q: Do a more general `listSplitSizes`?
    for trial in 1..3:              # ~10X overall slowdown factor w/defaults
      var rs: RunningStat
      for mi in subsampled(pnts, dimx, dimy, k, ni, score): rs.push mi
      varSum += F(ni - 1)/ni.F * rs.varianceS # Holmes2019Eq(9)seems to have..
      niSum += ni - 1                         #..missed canceling N from Eq(8).
  result[1] = sqrt(varSum/niSum.F)  #..Their Matlab _stddev code also drops it.

when isMainModule:
  import std/[syncio, strutils, enumerate], cligen, cligen/strUt
  when defined(danger): randomize()

  proc load(f: File, noise: F): seq[F] =
    var nCol = 0
    for i, row in enumerate(f.lines):
      let n0 = result.len
      for col in row.split:
        result.add col.strip.parseFloat + rand(1.0)*noise
      if nCol > 0 and result.len - n0 != nCol:
        stderr.write "stdin:",i,": irregular\n"
      nCol = result.len - n0

  proc estMI(k=6, file="", xDim=1, yDim=0, noise=0.0, score=false, err=4) =
    ##Estim.Mutual Info (bits) from data set of 1-line vector pairs by rectangle
    ##method of Kraskov, Stogbauer & Grassberger 2004, PhysRevE 69(6)066138 with
    ##follow-on ideas by Holmes & Nemenman 2019 PhysRevE 100(022404).
    let yDim = if yDim == 0: xDim else: yDim
    let pnts = load(if file.len > 0: file.open else: stdin, noise)
    let (val,err) = pnts.miErr(xDim, yDim, k, score); echo fmtUncertain(val,err)

  include cligen/mergeCfgEnv; dispatch estMI, help={
    "k": "k-th out nearest neighbor",
    "file": "filename; \"\" => stdin; Fmt: xDim+yDim rows",
    "xDim": "x vector dimension", "yDim": "y vector dimension (0=>same as x)",
    "noise": "U01-scale noise to add to each coord", "score": "t=normal score",
    "err": "max number of subdivisions to estim err with"}
