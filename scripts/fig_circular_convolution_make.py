#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure for Part XX: the circular convolution of castellated indicator functions
on T^1, and the three quantized interference states of Twin-Twin entanglement.
Requires: numpy, matplotlib.
"""
import math, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PRIMES = [5,7,11]
def holes(q): a=pow(6,-1,q); return sorted({a,(-a)%q})
def Rq(q,j):
    a=pow(6,-1,q); U={a,(-a)%q,(a-j)%q,(-a-j)%q}
    return (q-len(U))/q/((q-2)/q)**2

def draw_torus(ax, q, j, y, label):
    """Draw two castellated indicator bars (center N top, center N+j bottom) over [0,1)."""
    a = pow(6,-1,q); hA = {a,(-a)%q}; hB = {(a-j)%q,(-a-j)%q}
    # transmitting (white) vs blocked (dark) cells
    for r in range(q):
        # center N
        cA = "#d9d9d9" if r in hA else "#1f4e79"
        ax.add_patch(Rectangle((r/q, y+0.16), 1/q, 0.16, fc=cA, ec="w", lw=.6))
        # center N+j (shifted): blocked where r in hB
        cB = "#d9d9d9" if r in hB else "#b4341f"
        ax.add_patch(Rectangle((r/q, y), 1/q, 0.16, fc=cB, ec="w", lw=.6))
        # union-blocked marker (destructive overlap) -> hatch the transmitting cells that BOTH pass
    # transmitting overlap (both pass): green ticks under
    both = [r for r in range(q) if r not in hA and r not in hB]
    blocked = q - len(both)
    ax.text(1.02, y+0.08, f"{label}: blocked {blocked}/{q},  "
            r"$R_q=$"+f"{Rq(q,j):.3f}", va="center", fontsize=8.5)
    return blocked

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.0))
fig.suptitle("Part XX: circular convolution of castellated indicators, and the three interference states",
             fontsize=12.5, fontweight="bold")

# ---- Panel A: the three states on T^1 for q=7 ----
a = ax[0]; q = 7; inv3 = pow(3,-1,q)
states = [("State A  ($j\\equiv0$)", 0),
          ("State B  ($j\\equiv\\pm3^{-1}$)", inv3),
          ("State C  (otherwise)", 1)]
ys = [1.4, 0.7, 0.0]
for (lab, j), y in zip(states, ys):
    draw_torus(a, q, j, y, lab)
a.text(0.0, 2.02, r"center $N$ (blue) and shifted center $N+j$ (red);  grey = $\mathrm{dead}(q)$ hole",
       fontsize=8.5)
a.set_xlim(-0.02, 1.55); a.set_ylim(-0.15, 2.15)
a.set_xticks([0,0.5,1.0]); a.set_yticks([])
a.set_xlabel(r"phase $\theta\in[0,1)$  (residue / $q$),   $q=7$")
a.set_title("(A) two castellated bars on $\\mathbb{T}^1$; overlap = joint survival")

# ---- Panel B: R_q(j) quantization for q=5,7,11 ----
b = ax[1]; J = list(range(0,15))
col = {5:"#1f4e79",7:"#b4341f",11:"#2e7d32"}
for q in PRIMES:
    b.plot(J, [Rq(q,j) for j in J], "o-", c=col[q], ms=5, lw=1, label=f"$q={q}$")
b.axhline(1, color="k", ls="--", lw=.8)
b.text(14.2, 1.02, "indep.", fontsize=8, va="bottom", ha="right")
# annotate the three levels for q=5
for txt,val in [("A: q/(q-2)",5/3),("B",5*2/9),("C: q(q-4)/(q-2)²",5*1/9)]:
    b.annotate(txt,(0.1, val),fontsize=7.5,color="#1f4e79",va="center")
b.set_xlabel("lag $j$ (centre separation)"); b.set_ylabel(r"$R_q(j)$")
b.set_title("(B) three quantized states per prime (A/B $>1$, C $<1$)")
b.legend(fontsize=9); b.grid(alpha=.3)

fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig("fig_circular_convolution.png", dpi=200)
fig.savefig("fig_circular_convolution.pdf")
print("wrote fig_circular_convolution.{png,pdf}")
# print the global j=1 product (State C everywhere) for the caption
R1=1.0
for q in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]: R1*=Rq(q,1)
print(f"R(1)=prod R_q(1) (all State C) = {R1:.4f} = S_quad/S_twin^2")
