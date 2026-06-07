#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification for Part XX: the three quantized interference states of Twin-Twin
entanglement, and the spectral-cutoff k*(epsilon) error table.

For a prime q>3 the twin's two dead residues are {6^{-1}, -6^{-1}} (mod q). The
joint survival of twins at centres N and N+j is the circular-convolution overlap
of the two castellated indicators; for integer lag j it equals the discrete count
of co-surviving residues, giving the correlation factor

    R_q(j) = (q - |dead_A u dead_{B_j}|)/q / ((q-2)/q)^2 ,   dead_{B_j}=dead_A - j.

This script checks the closed three-state formula against that direct count, and
tabulates k*(epsilon): the smallest cutoff prime at which the truncated product
prod_{q<=k*} R_q(1) is within relative error epsilon of the full product.
Requires: nothing beyond the standard library.
"""
import math

def sieve_primes_gt3(limit):
    s = bytearray([1]) * (limit + 1); s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            for k in range(i*i, limit + 1, i): s[k] = 0
    return [p for p in range(5, limit + 1) if s[p]]

def dead_twin(q):
    a = pow(6, -1, q)
    return {a, (-a) % q}

def Rq_direct(q, j):
    """Direct union-overlap correlation factor."""
    A = dead_twin(q)
    Bj = {(r - j) % q for r in A}
    nu = len(A | Bj)
    return (q - nu) / q / ((q - 2) / q) ** 2

def Rq_formula(q, j):
    """Closed three-state formula."""
    inv3 = pow(3, -1, q)
    if j % q == 0:                              # State A: holes coincide
        return q / (q - 2)
    if j % q in (inv3 % q, (-inv3) % q):        # State B: one hole aligns
        return q * (q - 3) / (q - 2) ** 2
    return q * (q - 4) / (q - 2) ** 2           # State C: destructive

def verify_states():
    print("== three-state formula vs direct union count ==")
    ok = True
    for q in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59):
        for j in range(q):
            if abs(Rq_formula(q, j) - Rq_direct(q, j)) > 1e-12:
                ok = False
                print(f"  MISMATCH at q={q}, j={j}")
    print("  all states match:", ok)
    print("  sample (q=7): "
          + ", ".join(f"j={j}:{Rq_formula(7,j):.3f}" for j in range(7)))

def kstar_table(limit=200000):
    primes = sieve_primes_gt3(limit)
    ln_full = sum(math.log(Rq_formula(q, 1)) for q in primes)  # j=1: all State C
    print(f"\n== k*(epsilon) for the truncated product at j=1 ==")
    print(f"  full product R(1) = {math.exp(ln_full):.6f} = S_quad / S_twin^2")
    print(f"  {'epsilon':>9} {'k*':>8} {'primes kept':>12} {'error at k*':>13}")
    cum = 0.0; kept = 0
    targets = [1e-1, 1e-2, 1e-3, 1e-4]
    ti = 0
    rows = []
    for q in primes:
        cum += math.log(Rq_formula(q, 1)); kept += 1
        err = abs(math.exp(ln_full - cum) - 1)   # error of truncating AT q
        while ti < len(targets) and err < targets[ti]:
            rows.append((targets[ti], q, kept, err)); ti += 1
        if ti >= len(targets): break
    for eps, k, n, e in rows:
        print(f"  {eps:>9.0e} {k:>8} {n:>12} {e:>13.2e}")

if __name__ == "__main__":
    verify_states()
    kstar_table()
