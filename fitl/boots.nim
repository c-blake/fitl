##[ Define template to re-sample with replacement, aka statistical bootstrap.]##
import random; export random
if defined release: randomize()
type BootKind* = enum bMoving="moving", bCircular="circular" ## block boot mode

proc closestMultiple*(n, B: int): int =
  ## Return `n` rounded to the nearest multiple of `B`.  Stdlib addition?
  let L = (n div B)*B; let H = L + B
  if n - L < H - n: L else: H
#                                                  |-- Nim behavior blocks B=1
template forBoot*(n: int; boot=1000; mode=bMoving; B, put, fit) =
  ##[ Logic template to run `boot` successive estimations on data sets of size
  `n` re-sampled w/replacement in `mode bMoving|bCircular` blocks of length `B`.
  Caller provides `put(i, j)` to put 1 original point @`j` into some new data
  array@`i` (can be implicitly @end) & run a `fit`.  Centralizes `n mod B==0`
  check, `bMoving` end avoidance & `bCircular` wraparound.  While this only
  insists on n==K*B, it is statistically best for synthetic data set sizes to
  be as close to actual sample sizes as possible a la `closestMultiple`. ]##
  template putB(i, j0, B) =
    for j in j0 ..< j0 + B: put i, j

  if n mod B != 0: raise newException(ValueError,
    "`n` must be a multiple of block size `B`; Call with `closestMultiple`")
  for k in 1..boot:
    var i = 0                   # Re-sample a data set
    while i < n:
      let j0 = rand(n - (if mode == bMoving: B else: 1))
      case mode
      of bMoving: putB i, j0, B # Validity guaranteed by rand(n - B)
      of bCircular:             # Same as above, but with EOData wraparound
        if (let overflow = j0 + B - n; overflow > 0):
          putB i, j0, B - overflow
          putB i, 0 , overflow  # Copy B slots total: some above, rest here
        else:                   # Sample "kinda" needs S1 topology to make sense
          putB i, j0, B
      i += B
    fit                         # E.g., init aux, estimate, save params
