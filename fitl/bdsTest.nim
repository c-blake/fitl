## c-blake's translation of Blake LeBaron's C code to a Nim library.
import std/[strformat, algorithm, math], cligen/sysUt, cligen
template `+!`[T](p:ptr T, i:int):untyped=cast[ptr T](cast[int](p) +% i*T.sizeof)
template `-!`[T](p:ptr T, i:int):untyped=cast[ptr T](cast[int](p) -% i*T.sizeof)
template `-!`[T](q, p:ptr T):untyped=(cast[int](q) -% cast[int](p)) div T.sizeof
template mkItr(itr, cmp) =
  iterator itr[T](a, b: ptr T): ptr T =
    var i = a
    while cmp(cast[uint](i), cast[uint](b)): yield i; i = i +! 1
mkItr `..<`, `<`; mkItr `..`, `<=`

var   trace   = false
const NBIT    = 15 # Usable bits/word. 15 keeps lookup table small in count algo
const ALLBITS = 0xffff  #Q: Why is this not 0x7fff?
const TABLEN  = 32767
template bit(i): untyped = cast[int16](1 shl i)

proc gridOn(bStart: seq[ptr int16]; x, y: int) = # turn bits on
  if x != y:
    var ix = x; var iy = y
    if x > y: ix = y; iy = x
    iy = iy - ix - 1
    let ipos =       iy      div NBIT
    let ibit = NBIT - 1 - iy mod NBIT
    let p = bStart[ix] +! ipos
    p[] = p[] or ibit.bit

proc mkMask(r, n, nbit, drop: int; mask: pua int) =
  mask[1]=ALLBITS; mask[0]=mask[1]  # mask[0], mask[1]: 2-word mask.
  let last = (n - r - 1) div nbit   # Row `r`; `nbit`,`drop`: bits used&dropped
  for i in n - drop ..< n:
    let itrue = i - r - 1
    let j = last     - itrue div nbit
    let k = nbit - 1 - itrue mod nbit
    mask[j] = mask[j] xor k.bit

proc embed(bStart: seq[ptr int16]; n, dim: int) =
  for j in 0 ..< n - dim:     # Embed to next higher dim; `g(i,j) &= g(i+1,j+1)`
    var i = bStart[j]
    for i2 in bStart[j + 1] ..< bStart[j + 2]:
      i[] = i[] and i2[]; i = i +! 1
    if i != bStart[j + 1]: i[] = 0

type
  Pos[F] = tuple[val: F, pos: int]
  BDS*[F] = object      ## Holds state for BDS stat calc for `n` observations.
    posTab: seq[Pos[F]]
    mask, lookup: seq[int]
    grid  : seq[int16]
    start : seq[ptr int16]

proc space[T](x: seq[T]): int = x.len*x[0].sizeof
proc space*[F](b: BDS[F]): int = ## Lower bound indirect space usage of `b`
  b.sizeof+b.posTab.space+b.mask.space+b.lookup.space+b.grid.space+b.start.space

proc initBDS*[F](n: int): BDS[F] = ## Make BDS calc holder for `n` observations.
  result.posTab.setLen n
  result.mask.setLen   2*n
  result.lookup.setLen TABLEN + 1
  result.start.setLen  n + 1
  var sz = 0                            # Find grid size
  for i in 0..n: sz.inc (n - i) div NBIT  +  1
  result.grid.setLen sz                 # Grid is defined as short (2 byte ints)
  result.start[0] = result.grid[0].addr
  for i in 1..n: result.start[i] = result.start[i - 1] +! ((n - i) div NBIT + 1)
  for i in 0..TABLEN:                   # Table for bit counting
    for j in 0 ..< NBIT: (if (i and j.bit) != 0: inc result.lookup[i])

proc evalc[F](b: BDS[F]; n: int): float =
  var cnt = 0   # Return count stats for grid; Zero uncounted parts using mask
  for j in 0..<n:
    if b.start[j + 1] -! b.start[j] > 2:
      for i in b.start[j] ..< b.start[j + 1] -! 2:
        cnt.inc b.lookup[i[]]
        if b.lookup[i[]] > NBIT: ResourceExhausted !! &"{i[]} {b.lookup[i[]]}"
      for i in b.start[j + 1] -! 2 ..< b.start[j + 1]:
        cnt.inc b.lookup[i[] and b.mask[2*j + (b.start[j + 1] -! i - 1)]]
    else:
      for i in b.start[j] ..< b.start[j + 1]:
        cnt.inc b.lookup[i[] and b.mask[2*j + (b.start[j + 1] -! i - 1)]]
  if trace: echo "count = ",cnt
  return 2*cnt.float / float(n*(n - 1))

proc kc*[F](b: var BDS[F]; x: openArray[F]; m=2; drop = -1; eps: F):
       tuple[k: float, c: seq[float]] =
  ## `x`: time series to test;  Length must be same as in `b=initBDS(x.len)`.
  ## `m`: cstats will be done for dim/lag 1..m;  `drop`: num.data to drop @end.
  ##   { c(2) can use more data than c(3),.. => Ease ignoring last few so all
  ##     c() are done over same data.  Eg. for m=3 we might use x(1..n-2). }
  ## `eps`: epsilon value for close points.
  ## Returns `k` & `c`: 1-origin indexed raw c values c[1], c[2], c[3], ...
  if x.len != b.posTab.len: Value !! &"{x.len} != {b.posTab.len}"
  let drop = if drop != -1: drop else: m - 1  # -1=>All c[] estim w/n-m+1 points
  let n = x.len; let nOb = n - drop; let nObF = nOb.float
  for ip in b.grid[0].addr ..< b.start[n] +! 1: ip[] = 0
  for i in 0..<n: b.posTab[i] = (x[i], i)
  b.posTab.sort
  let posTLast = b.posTab[n - 1].addr
  var count = 0                         # Row-by-row build; Use Theiler method
  var phi = 0.0
  for p in b.posTab[0].addr .. posTLast:
    var tcount = 0
    var pt = p
    while pt.val - p.val <= eps:        # Count to right
      b.start.gridOn p.pos, pt.pos
      if p.pos < nOb and pt.pos < nOb: inc tcount
      if pt == posTLast: break
      else: pt = pt +! 1
    if p != b.posTab[0].addr:           # Count to left; For Dechert k NOT grid
      pt = p -! 1
      while p.val - pt.val <= eps:
        if p.pos < nOb and pt.pos < nOb: inc tcount
        if pt == b.posTab[0].addr: break
        else: pt = pt -! 1
    count.inc tcount
    phi += float(tcount*tcount)         # Dechert speed up K
  count -= nOb                          # Adjust k and c to u statistic
  phi   -= float(nOb + 3*count)
  if trace: echo &"{count} {phi}"
  result.k = phi/(nObF*(nObF - 1)*(nObF - 2))
  result.c.setLen m + 1
  result.c[1] = count.float / (nObF*(nObF - 1))
  for i in 0..<nOb: mkMask(i, n, NBIT, drop, b.mask[2*i].addr.toPua) # Make Mask
  for i in 2..m: b.start.embed n, i; result.c[i] = b.evalc nOb       # HigherDim

proc cstat*(c, cm, k: float; m, n: int): float =
  ## c&k->test stat asymptotically~N(0,1).  Brock,Hsieh,LeBaron 1991 Chap2,pg43.
  ## `c`: c[1]; `cm`: c[m]; `k`: k stat; `m`: embedding dim; `n` = nPoints.
  var v = 0.0; let M = m.float
  for j in 1..<m: v += 2*k^(m - j)*c^(2*j)
  v += k^m + (M - 1)*(M - 1)*c^(2*m) - M*M*k*c^(2*m - 2)
  (cm - c^m)*sqrt(n.float/v)*0.5

proc sampleVar(x: openArray[float]): float =
  var s, s2: float; let n = x.len.float
  for e in x: s += e; s2 += e*e
  (s2/n - (s/n)*(s/n))*n/(n - 1.0)

proc BDS_N01s*(x: openArray[float]; m: int; es=0.5): seq[float] =
  ## A simple entry point, much like the R `tseries::bdstest`. Returns No.Sigmas
  ## for each lag tuple size 2..`m`.  Asymptopia at n >~ 500.
  var b = initBDS[float](x.len)         # Raw k&c statistics
  let (k, c) = b.kc(x, m, m - 1, es*sqrt(x.sampleVar))
  if trace: (echo(&"k = ",k); for i in 1..m: echo &"c({i}) {c[i]}")
  result.setLen m + 1
  for i in 2..m: result[i] = cstat(c[1], c[i], k, i, x.len - m + 1)
  if trace: echo "b.space: ",b.space," bytes"

when isMainModule:
  proc test(x: seq[float]; m=3; es=0.5; trace=false) =
    ## BDS Test command (like R tseries::bdstest).  Asymptopia at n >~ 500.
    if x.len < m + 1: Value !! "\e[1mNOT ENOUGH DATA\e[m"
    bdsTest.trace = trace; for e in BDS_N01s(x, m, es)[2..^1]: echo &"{e:.2f}"
  dispatch test, cmdName="bdsTest"      # 1.1 2.05 3.03 4.02 5.01 6.0 7.1 8.05
