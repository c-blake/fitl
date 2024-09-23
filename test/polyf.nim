# E.g.: nim r test/gen -g0.1 | awk '{print $3,$1}' | nim r test/polyf
when not declared(stdin): import std/[syncio, formatfloat]
import fitl/svdx, std/strutils

var x, y, w, b: seq[float]
for line in stdin.lines:
  let cols = line.split()
  x.add parseFloat(cols[0])
  y.add parseFloat(cols[1])
  if cols.len > 2: w.add parseFloat(cols[2])

b.polyFit x, y, w
echo b
