import math, basicLA
{.passC: "-O3 -ffast-math -march=native -mtune=native".}

proc mean[F](x: ptr F; n: int): F = sum(x, n) / F(n)

proc corr[F](x, y: openArray[F]; xm, ym: F; n: int): F =
  var syy, sxy, sxx: F
  for i in 0 ..< n:
    sxx += (x[i] - xm)*(x[i] - xm)
    syy += (y[i] - ym)*(y[i] - ym)
    sxy += (x[i] - xm)*(y[i] - ym)
  sxy / sqrt(sxx * syy)

proc corr*[F](x, y: openArray[F]): F =
  ## Return Pearson linear correlation coefficient between x & y (over leading
  ## elements if one array is longer).
  let n = min(x.len, y.len)
  corr x, y, mean(x[0].unsafeAddr, n), mean(y[0].unsafeAddr, n), n

proc corrAuto*[F](x: openArray[F], lag=1): F =
  if lag < x.len: corr(x, x[lag..^1]) else: F(0)

proc covMat*[F](v: var openArray[F]; x: openArray[F]; n, m: int) =
  ##Save in v usual symmetric m*m Covariance matrix for n*m input matrix x where
  ##x[i+n*j] is the j-th column of the i-th row.  I.e., samples is the faster
  ##moving index, not variables aka column-major.
  var ni = if n > 1: 1.0 / F(n-1) else: 1.0
  var xm = newSeq[F](m) # means
  for j in 0 ..< m: xm[j] = mean(x[n*j].unsafeAddr, n)
  for i in 0 ..< m:     # lower triangle & diag
    for j in 0..i:
      v[m*i+j] = dots(x[n*i].unsafeAddr, x[n*j].unsafeAddr, xm[i], xm[j], n)*ni
  for i in 0 ..< m:     # fill in upper from lower
    for j in i + 1 ..< m: v[m*i+j] = v[m*j+i]

proc covMat*[F](x: openArray[F]; m, n: int): seq[F] =
  ## Return usual symmetric m*m Covariance matrix for n*m input matrix x where
  ## x[i+n*j] is the j-th column of the i-th row.  I.e., samples is the faster
  ## moving index, not variables.
  result.setLen m*m
  covMat(result, x, m, n)

when isMainModule:
  echo corrAuto([1.0, 1.9, 3.2, 4.1])
  echo covMat([ 1.0, 2, 3, 4, 5, 6, 7, 8 ], 4, 2)
