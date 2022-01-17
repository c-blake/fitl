## This is just a few BLAS Level 1 operations.  With gcc -ffast-math auto-vec
## they are about as fast (or sometimes a bit faster) than optimized libs. (A
## bit faster since they are more specific to dense/non-strided iterations.)
{.passc: "-O3 -ffast-math -march=native -mtune=native".}

proc sum*[F](x: ptr F; n: int): F =
  ## Sum of elements with accumulator in same arithmetic width as params.
  let x = cast[ptr UncheckedArray[F]](x)
  for i in 0..<n: result += x[i]

template sum0*(i, lim, ex: untyped): untyped =
  ## Sum of expressions ex(i) from 0..<`lim`.  Accumulator is `type(ex)`.  Eg.,
  ## `sum0(i, 5, float64(i))==10.0` accumulated in float64 arithmetic.
  block:
    var i {.inject.}: type(lim)   # put i in scope for ex(i); init to 0
    var tot: type(ex)             # accumulator type taken from ex
    for i in i ..< lim: tot += ex # make i same type as lim & read-only within ex
    tot

template Sum*[T](i; itr: iterable[T]; ex): untyped =
  ## A more general version of `sum0` leveraging Nim iterator builder syntax to
  ## mimick math notation.  Eg., `Sum(i, 1..10, 2*i)`.  It is capitalized to be
  ## like math & to not collide w/math.sum (untyped BREAKS OVERLOAD RESOLUTION).
  block:
    var i {.inject.}: T # put i in scope for ex(i)
    var tot: type(ex)   # accumulator type from ex(i)
    for kji in itr:     # engage generic iterator yielding T
      let i = kji       # rename itr yield to read-only i
      tot += ex
    tot

proc dot*[F](x, y: ptr F; n: int): F {.inline.} =
  ## Dot product with accumulator in same arithmetic width as params.
  let x = cast[ptr UncheckedArray[F]](x)
  let y = cast[ptr UncheckedArray[F]](y)
  for i in 0 ..< n: result += x[i] * y[i]

proc dots*[F](x, y: ptr F; x0, y0: F; n: int): F {.inline.} =
  ## Dot (x-x0).(y-y0) with accumulator in same arithmetic width as params.
  let x = cast[ptr UncheckedArray[F]](x)
  let y = cast[ptr UncheckedArray[F]](y)
  for i in 0 ..< n: result += (x[i] - x0) * (y[i] - y0)

proc xpose*[F](t: var openArray[F]; x: openArray[F]; n, m: int) =
  ## Fill pre-sized `t` with the transpose of n*m matrix `x`.
  let t = cast[ptr UncheckedArray[F]](t)
  let x = cast[ptr UncheckedArray[F]](x)
  for j in 0 ..< m: # Order for dense wr sinces rd has ok chance of being cached
    for i in 0 ..< n: t[n*j + i] = x[m*i + j]

proc xpose*[F](x: openArray[F]; n, m: int): seq[F] =
  ## Return seq-allocated transpose of n*m matrix `x`.
  result.setLen n*m
  xpose(result, x, n, m)

proc vadd*[F](x: ptr F; n: int, v: F) {.inline.} =
  ## Add `v` to each element of `x`
  let x = cast[ptr UncheckedArray[F]](x)
  for i in 0 ..< n: x[i] += v

proc vmul*[F](x: ptr F; n: int, v: F) {.inline.} =
  ## Multiply each element of `x` by `v`.
  let x = cast[ptr UncheckedArray[F]](x)
  for i in 0 ..< n: x[i] *= v

func mnmx*[F: SomeFloat](xs: openArray[F]): (F, F) =
  ## One pass data range (min, max).
  if xs.len > 0:
    result[0] = xs[0]; result[1] = xs[0]
    for x in xs: result[0] = min(result[0], x); result[1] = max(result[1], x)
  else: result[0]=F(NaN); result[1]=F(NaN)

func mvar*[F: SomeFloat](xs: openArray[F]): (F, F) =
  ## One pass mean & variance; NOTE: pop variance (1/n).
  if xs.len > 0: # --passC:-ffast-math can autovec whole loop into vsubp[sd]/vaddp[sd] insns
    let dx = xs[0]
    var av, vr: F
    for x in xs:
      let x = x - dx; av += x; vr += x*x
    av /= F(xs.len); vr /= F(xs.len)
    result[0] = F(av.float + dx)
    result[1] = F(max(1e-30*av*av, vr - av*av))
  else: result[0]=F(NaN); result[1]=F(NaN)

func mvars*[F: SomeFloat](xs: openArray[F]): (F, F) =
  ## One pass mean & sample variance
  let mV = xs.mvar
  (mV[0], mV[1]*F(xs.len)/F(xs.len - 1))

proc newSeq*[T](len: Natural, val: T): seq[T] =
  ## Allocate & initialize a seq to a value
  result.setLen len
  for i in 0..<len: result[i] = val

proc zero*[T](x: var openArray[T]) =
  ## Zero memory of `x`.
  zeroMem x[0].addr, x.len*T.sizeof

proc set*[T](x: var openArray[T], v: T) =
  ## Set every `x[i]` to `v`.
  for e in mitems(x): e = v
