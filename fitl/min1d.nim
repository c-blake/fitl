import math # Module with algos to find minima of 1-D scalar functions.

type ScalarF*[F] = proc(x: F): F
  ## Nim requires casting proc types too much.

proc minBrent*[F](x: var F; a, b: F; f: ScalarF[F]; tol=1e-7, itMx=100,
                  dxMin0=1e-9): F {.discardable.} = # can want only arg min
  ## Return minimum of a univariate function f(x) on [a,b] by Brent's method.
  ## x is also set to arg min_[a,b] f(x).
  const golden = F(0.38196601125010515180)  # 1-1/phi where phi^2-phi==1.
  var a = a; var b = b
  var u, v, w,  fx, fu, fv, fw,  xm, dx0, dxMin, dx, s: F
  x  = a + golden*(b - a); w  = x ; v  = x  # Initial evaluation
  fx = f(x)              ; fw = fx; fv = fx
  for it in 1..itMx:
    xm    = (a + b)/F(2)
    dxMin = tol*abs(x) + F(dxMin0)
    if abs(x - xm) < F(2)*dxMin - (b - a)/F(2):
      break                                 # Converged
    if abs(dx) <= dxMin:                    # Golden section step
      dx = (if x >= xm: a else: b) - x
      s  = golden*dx
    else:                                   # Parabolic interpolation step
      let xw = (x - w)*(fx - fv)
      let xv = (x - v)*(fx - fw)
      var vw = F(2)*(xv - xw)
      var p  = (x - v)*xv - (x - w)*xw
      if vw > F(0): p  = -p
      else        : vw = -vw
      dx0 = dx
      dx  = s
      if p > vw*(a - x) and p < vw*(b - x) and abs(p) < abs(vw*dx0/F(2)):
        s = p / vw                          # Parabolic step is useful
        u = x + s
        if (u - a) < F(2)*dxMin or (b - u) < F(2)*dxMin:
          s = copySign(dxMin, xm - x)
      else:                                 # Unuseful => golden section
        dx = (if x >= xm: a else: b) - x
        s  = golden*dx
    u  = x + (if abs(s) >= dxMin: s else: copySign(dxMin, s))
    fu = f(u)                               # Split by at least dxMin
    if fu > fx:
      if u < x: a = u
      else    : b = u
      if fu <= fw or w == x:
        v = w; fv = fw
        w = u; fw = fu
      elif fu <= fv or v == x or v == w:
        v = u; fv = fu
    else:
      if u >= x: a = x
      else     : b = x
      v = w; fv = fw
      w = x; fw = fx
      x = u; fx = fu
  fx

when isMainModule:
  when not declared(addFloat): import std/formatfloat
  var nEval = 0
  proc g(x: float): float = (inc nEval; ln(1.0 + (x - 2.0)*(x - 2.0)))
  var x: float
  echo "min: ", minBrent(x, 0.0, 100.0, g), " at: ",x," after ",nEval," evals"
