## k-sample statistical homogeneity test based on Fritz Scholz' R package.
{.warning[Uninit]:off, warning[ProveInit]:off.} # Should be verbosity:2 not 1
import std/[math, random, algorithm]
type F = float
func count(dat: openArray[F], z: F): int = (for e in dat: (if e==z: inc result))
func getSum(x: openArray[int], lim=0): int = (for i in 0..<lim: result += x[i])
func dedup(a: var seq[F]) =
  var i = 0
  for j in 1 ..< a.len:
    if a[j] != a[i]: inc i; a[i] = a[j]
  a.setLen i + 1

func sortU*(x: var seq[F]) = x.sort; x.dedup

func adkStat*(x, u: openArray[F]; ns: openArray[int]): (F, F) =
  ## Return (orig, alt) k-sample Anderson-Darling test stats for non-parametric
  ## (rank) test of Scholz & Stephens 1987.  Here `x` is the concatenated sample
  ## while `ns` are subsample sizes in the same order.  `u` holds the L distinct
  ## ordered observations in the pooled sample.
  let (k, L) = (ns.len, u.len)  # Q: if L == 1: return (0.0, 0.0)
  var fij  = newSeq[int](k*L)   # num.pts in ith samp for u[j] where i=1-k,j=1-L
  var lvec = newSeq[int](L)     # lvec[j] = multiplicity of u[j]
  var samp = newSeq[seq[F]](k)  # Samples being compared
  var N = 0
  for i, n in ns:
    samp[i].setLen n
    for j in 0..<n: samp[i][j] = x[N + j]
    N += ns[i]
  for j in 0..<L:
    for i in 0..<k:
      fij[i + j*k] = samp[i].count(u[j])
      lvec[j] += fij[i + j*k]
  for i in 0..<k:
    var mij, maij, innSum, aInnSum: F
    for j in 0..<L:
      mij += fij[i + j*k].F
      maij = mij - 0.5*fij[i + j*k].F
      let bj = getSum(lvec, j + 1).F
      let baj = bj - 0.5*lvec[j].F
      if j < L - 1:
        let t = N.F*mij - ns[i].F*bj
        innSum += lvec[j].F*t*t/(bj*(N.F - bj))
      let t = N.F*maij - baj*ns[i].F
      aInnSum += lvec[j].F*t*t/(baj*(N.F - baj) - N.F*lvec[j].F*0.25)
    result[0] += innSum/ns[i].F         # AkN2
    result[1] += aInnSum/ns[i].F        # AakN2
  let nInv = 1.0/N.F                    # (orig,alt) k-sample A-D test stats
  result[0] = result[0]*nInv            # AkN2
  result[1] = (1 - nInv)*result[1]*nInv # AakN2

const N: ptr seq[F] = nil
proc adkSim*(x: openArray[F]; ns: openArray[int]; m=99; s0=N, sA=N): (int, int)=
  ## Count by-chance (orig,alt) k-sample Anderson-Darling test stats >= observed
  ## for non-parametric (rank) test described in Scholz,Stephens 1987 for later
  ## `binomP` p-value estimates.  Here, `x` is a catenated list of observations,
  ## `ns` are the k sample sizes (in the same order), `m` is the number of sims,
  ## and `s0, sA` optional saved test statistic values.
  var u = x[0..^1]                      # Copy x for sort & de-duplicate
  var x = u                             # Copy x again for shuffling
  u.sortU                               # sort & de-duplicate
  let (datA0, datAA) = x.adkStat(u, ns) # Get observed AkN2,AakN2
# echo "A0: ",datA0," AA: ",datAA
  for i in 1..m:                        # Simulate
    x.shuffle
    let (a1, aa) = x.adkStat(u, ns)     # Get simulated AkN2,AakN2
    if not s0.isNil: s0[].add a1
    if not sA.isNil: sA[].add aa
    inc result[0], (a1 >= datA0).int
    inc result[1], (aa >= datAA).int

when isMainModule:
 when defined danger: randomize()
 let x=[0.824,0.216,0.538,0.685,0.448,0.348,0.443,0.722,0.403,0.268,0.440,0.087]
 echo adkSim(x, [4,4,4], m=1000000) # ADStatVsn2 = 2.6240 != 2.6200, pV matches
