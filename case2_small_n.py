#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:58:13 2026

@author: chiarabellotti
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from fractions import Fraction
from math import isqrt
import sys


# ============================================================
# Settings
# ============================================================

N_START = 50
N_END = 199
TARGET = Fraction(1, 3)


# ============================================================
# Basic arithmetic
# ============================================================

def primes_up_to(limit: int):
    """Return all primes <= limit."""
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [p for p in range(2, limit + 1) if sieve[p]]


def Q_n(n: int, primes):
    """
    Q_n = {p prime : n/40 < p <= n-3}.

    We use 40*p > n instead of p > n/40 to avoid floating-point issues.
    """
    return [p for p in primes if 40 * p > n and p <= n - 3]


# ============================================================
# Exact polynomial arithmetic
# ============================================================

def factor_poly(step: int, alpha: Fraction, degree: int):
    """
    Return the polynomial coefficients up to x^degree for

        (1 - x^step)^alpha.

    The result is a list coeffs such that

        coeffs[i] = coefficient of x^i.

    Coefficients are exact Fractions.
    """
    coeffs = [Fraction(0, 1) for _ in range(degree + 1)]

    # Coefficient of y^j in (1-y)^alpha is binom(alpha,j)(-1)^j.
    c = Fraction(1, 1)
    j = 0

    while j * step <= degree:
        coeffs[j * step] = c

        j += 1
        # Update:
        # c_j = c_{j-1} * (alpha - (j-1)) / j * (-1)
        c *= Fraction(alpha - (j - 1), j) * (-1)

    return coeffs


def multiply_truncated(A, B, degree: int):
    """
    Multiply two polynomials, keeping only terms up to x^degree.
    Coefficients are Fractions.
    """
    C = [Fraction(0, 1) for _ in range(degree + 1)]

    for i, ai in enumerate(A):
        if ai == 0:
            continue

        max_j = degree - i
        for j in range(max_j + 1):
            bj = B[j]
            if bj != 0:
                C[i + j] += ai * bj

    return C


def coeff_with_one_over_one_minus_x(factors, degree: int):
    """
    Compute exactly

        [x^degree] product(factors) / (1-x),

    where each factor is of the form

        (1 - x^step)^alpha.

    Input:
        factors = [(step_1, alpha_1), ..., (step_r, alpha_r)].

    Since multiplying by 1/(1-x) takes partial sums, if P(x) is the product
    of the factors, then

        [x^degree] P(x)/(1-x) = sum_{i=0}^degree [x^i]P(x).
    """
    poly = [Fraction(0, 1) for _ in range(degree + 1)]
    poly[0] = Fraction(1, 1)

    for step, alpha in factors:
        F = factor_poly(step, alpha, degree)
        poly = multiply_truncated(poly, F, degree)

    return sum(poly[: degree + 1])


# ============================================================
# Exact probabilities
# ============================================================

def prob_A_p(n: int, p: int):
    """
    Exact value of

        Pr(A_p)
        =
        1/p * [x^{n-p}] (1-x^p)^{1/p}/(1-x).
    """
    degree = n - p

    coeff = coeff_with_one_over_one_minus_x(
        factors=[(p, Fraction(1, p))],
        degree=degree,
    )

    return Fraction(1, p) * coeff


def prob_A_p_inter_A_q(n: int, p: int, q: int):
    """
    Exact value of

        Pr(A_p cap A_q)
        =
        1/(pq) * [x^{n-p-q}]
        (1-x^p)^{1/p}
        (1-x^q)^{1/q}
        (1-x^{pq})^{-1/(pq)}
        /(1-x).

    If p+q > n, the probability is zero.
    """
    if p + q > n:
        return Fraction(0, 1)

    degree = n - p - q

    coeff = coeff_with_one_over_one_minus_x(
        factors=[
            (p, Fraction(1, p)),
            (q, Fraction(1, q)),
            (p * q, Fraction(-1, p * q)),
        ],
        degree=degree,
    )

    return Fraction(1, p * q) * coeff


# ============================================================
# Exact second moment bound
# ============================================================

def exact_second_moment_bound(n: int, primes):
    """
    Compute exactly

        (M_n^*)^2 / (M_n^* + 2 R_n^*),

    where

        M_n^* = sum_p Pr(A_p)

    and

        R_n^* = sum_{p<q} Pr(A_p cap A_q).

    Here p,q run over Q_n = {p prime : n/40 < p <= n-3}.
    """
    Q = Q_n(n, primes)

    M_star = Fraction(0, 1)
    for p in Q:
        M_star += prob_A_p(n, p)

    R_star = Fraction(0, 1)
    for i, p in enumerate(Q):
        for q in Q[i + 1 :]:
            if p + q <= n:
                R_star += prob_A_p_inter_A_q(n, p, q)

    if M_star == 0:
        return Fraction(0, 1), M_star, R_star, len(Q)

    value = (M_star * M_star) / (M_star + 2 * R_star)

    return value, M_star, R_star, len(Q)


# ============================================================
# Main verification
# ============================================================

def verify_range():
    primes = primes_up_to(N_END)

    print(f"Checking exact second-moment bound for {N_START} <= n <= {N_END}")
    print(f"Target = {TARGET}")
    print()

    min_value = None
    min_n = None
    min_M = None
    min_R = None
    min_Q_size = None

    for n in range(N_START, N_END + 1):
        value, M_star, R_star, q_size = exact_second_moment_bound(n, primes)

        if min_value is None or value < min_value:
            min_value = value
            min_n = n
            min_M = M_star
            min_R = R_star
            min_Q_size = q_size

        if value <= TARGET:
            print("Verification failed.")
            print(f"First failing n = {n}")
            print(f"lower bound exact = {value}")
            print(f"lower bound decimal = {float(value):.15f}")
            print(f"M_star = {M_star}")
            print(f"R_star = {R_star}")
            print(f"|Q_n| = {q_size}")
            sys.exit(1)

    print("Verification successful.")
    print(f"For every {N_START} <= n <= {N_END},")
    print("the exact second-moment lower bound is > 1/3.")
    print()
    print(f"Minimum attained at n = {min_n}")
    print(f"Minimum exact value = {min_value}")
    print(f"Minimum decimal value = {float(min_value):.15f}")
    print(f"M_star at minimum = {min_M}")
    print(f"R_star at minimum = {min_R}")
    print(f"|Q_n| at minimum = {min_Q_size}")


if __name__ == "__main__":
    verify_range()