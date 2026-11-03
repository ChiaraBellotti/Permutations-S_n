#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:51:54 2026

@author: chiarabellotti
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import sys
from bisect import bisect_left, bisect_right

# ============================================================
# Settings
# ============================================================

N_START = 200
N_END = 9052
TARGET = 1.0 / 3.0

# We use a small tolerance only to avoid floating-point edge issues.
# The actual minimum is far above 1/3.
EPS = 1e-14


# ============================================================
# Basic arithmetic
# ============================================================

def primes_up_to(limit: int):
    """Return all primes <= limit."""
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [p for p in range(2, limit + 1) if sieve[p]]


def harmonic_numbers(limit: int):
    """
    Return list H with H[m] = 1 + 1/2 + ... + 1/m,
    and H[0] = 0.
    """
    H = [0.0] * (limit + 1)
    s = 0.0

    for m in range(1, limit + 1):
        s += 1.0 / m
        H[m] = s

    return H


def reciprocal_prime_prefix(primes):
    """
    prefix[i] = sum_{j < i} 1 / primes[j].
    """
    prefix = [0.0]
    s = 0.0

    for p in primes:
        s += 1.0 / p
        prefix.append(s)

    return prefix


# ============================================================
# M_n and R_n
# ============================================================

def primes_in_Q(n: int, primes):
    """
    Q_n = {p prime : n/40 < p <= n-3}.
    """
    lo = n / 40.0
    hi = n - 3

    i = bisect_right(primes, lo)
    j = bisect_right(primes, hi)

    return primes[i:j], i, j


def compute_M(n: int, Q, H):
    """
    M_n =
      sum_{n/40 < p <= n-3}
        1/p * (1 - H_{floor((n-p)/p)} / p).
    """
    total = 0.0

    for p in Q:
        r = (n - p) // p
        total += (1.0 / p) * (1.0 - H[r] / p)

    return total


def compute_R(n: int, Q):
    """
    R_n =
      sum_{p<q, n/40<p<q<=n-3, p+q<=n} 1/(pq).

    This is computed in O(|Q| log |Q|), using prefix sums of 1/p.
    """
    if not Q:
        return 0.0

    prefix = [0.0]
    s = 0.0

    for p in Q:
        s += 1.0 / p
        prefix.append(s)

    total = 0.0

    for j, q in enumerate(Q):
        # We need p < q and p <= n-q.
        max_p = min(q - 1, n - q)

        # Number of primes p in Q with p <= max_p.
        k = bisect_right(Q, max_p)

        # Since q is Q[j], p<q means p occurs before index j.
        k = min(k, j)

        if k > 0:
            sum_recip_p = prefix[k]
            total += (1.0 / q) * sum_recip_p

    return total


def lower_bound(n: int, primes, H):
    """
    Compute
        M_n^2 / (M_n + 2 R_n).
    """
    Q, _, _ = primes_in_Q(n, primes)

    M = compute_M(n, Q, H)
    R = compute_R(n, Q)

    if M <= 0:
        return 0.0, M, R, len(Q)

    value = (M * M) / (M + 2.0 * R)

    return value, M, R, len(Q)


# ============================================================
# Main verification
# ============================================================

def verify_range():
    primes = primes_up_to(N_END)
    H = harmonic_numbers(N_END)

    print(f"Checking {N_START} <= n <= {N_END}")
    print(f"Target = {TARGET:.15f}")
    print()

    best_value = float("inf")
    best_n = None
    best_M = None
    best_R = None
    best_q_size = None

    for n in range(N_START, N_END + 1):
        value, M, R, q_size = lower_bound(n, primes, H)

        if value < best_value:
            best_value = value
            best_n = n
            best_M = M
            best_R = R
            best_q_size = q_size

        if value <= TARGET + EPS:
            print("Verification failed.")
            print(f"First failing n = {n}")
            print(f"Lower bound = {value:.15f}")
            print(f"M_n = {M:.15f}")
            print(f"R_n = {R:.15f}")
            print(f"|Q_n| = {q_size}")
            sys.exit(1)

    print("Verification successful.")
    print(f"For every {N_START} <= n <= {N_END},")
    print("M_n^2 / (M_n + 2 R_n) > 1/3.")
    print()
    print(f"Minimum value = {best_value:.15f}")
    print(f"Attained at n = {best_n}")
    print(f"M_n = {best_M:.15f}")
    print(f"R_n = {best_R:.15f}")
    print(f"|Q_n| = {best_q_size}")


if __name__ == "__main__":
    verify_range()