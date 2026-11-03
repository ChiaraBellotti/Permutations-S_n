#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:42:27 2026

@author: chiarabellotti
"""

from collections import Counter
from fractions import Fraction
from math import factorial, isqrt
from typing import Optional


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


def partitions_of_n(n: int, max_part: Optional[int] = None):
    """
    Generate all integer partitions of n as tuples in weakly decreasing order.

    Example:
        partitions_of_n(4) yields
        (4,), (3,1), (2,2), (2,1,1), (1,1,1,1)
    """
    if max_part is None or max_part > n:
        max_part = n

    if n == 0:
        yield ()
        return

    for first in range(max_part, 0, -1):
        for rest in partitions_of_n(n - first, first):
            yield (first,) + rest


def class_weight(partition):
    """
    Exact probability weight of a conjugacy class with given cycle partition.

    The weight is
        1 / C(lambda),
    where
        C(lambda) = prod_k k^{m_k} m_k!.
    """
    counts = Counter(partition)

    denom = 1
    for k, m in counts.items():
        denom *= (k ** m) * factorial(m)

    return Fraction(1, denom)


def is_pre_p_cycle_partition(partition, p: int) -> bool:
    """
    For prime p, a partition corresponds to a pre-p-cycle iff

        m_p = 1
        m_{ip} = 0 for all i >= 2.
    """
    counts = Counter(partition)

    if counts.get(p, 0) != 1:
        return False

    total_n = sum(partition)

    multiple = 2 * p
    while multiple <= total_n:
        if counts.get(multiple, 0) != 0:
            return False
        multiple += p

    return True


def is_pre_prime_cycle_partition(partition, primes) -> bool:
    """
    True iff the partition is a pre-p-cycle for at least one prime p
    with 2 <= p <= n-3.
    """
    for p in primes:
        if is_pre_p_cycle_partition(partition, p):
            return True

    return False


def pi_n(n: int):
    """
    Exact value of pi_n: the proportion of elements of S_n that are
    pre-p-cycles for some prime p with 2 <= p <= n-3.
    """
    if n < 5:
        return Fraction(0, 1)

    primes = primes_up_to(n - 3)
    total = Fraction(0, 1)

    for part in partitions_of_n(n):
        if is_pre_prime_cycle_partition(part, primes):
            total += class_weight(part)

    return total


if __name__ == "__main__":
    target = Fraction(1, 3)

    print(" n   pi_n exact                         pi_n decimal")
    print("---------------------------------------------------------")

    all_ok = True

    for n in range(5, 50):
        value = pi_n(n)

        print(f"{n:2d}   {str(value):30s}   {float(value):.15f}")

        if 8 <= n <= 49 and value <= target:
            all_ok = False

            print()
            print("Verification failed:")
            print(f"n = {n}")
            print(f"pi_n = {value} = {float(value):.15f}")
            print(f"1/3  = {target} = {float(target):.15f}")

            break

    print()

    if all_ok:
        print("Verification successful:")
        print("pi_n >= 1/3 for all 8 <= n <= 49.")