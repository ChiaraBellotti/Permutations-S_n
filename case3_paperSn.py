"""Verify Case 3
"""

from array import array
from fractions import Fraction
from math import exp, isqrt
import sys

N_START = 9053
N_END = 1_199_999
A = 12
Q = 10**18  # Fixed-point scaling factor
TARGET_NUM = 333777
TARGET_DEN = 1_000_000
PROGRESS_EVERY = 100_000

# This strict rational upper bound on exp(Euler's constant) is proved below.
EGAMMA_NUM = 1_781_072_418
EGAMMA_DEN = 1_000_000_000
EULER_GAMMA_APPROX = 0.5772156649015328606  # PRINTING ONLY


def ceil_div(numerator, denominator):
    return (numerator + denominator - 1) // denominator


def check_euler_constant_bound():
    N = 2048
    harmonic_upper = sum(ceil_div(Q, j) for j in range(1, N + 1))
    log2_lower = sum(
        (2 * Q) // ((2 * j + 1) * 3 ** (2 * j + 1))
        for j in range(25)
    )
    gamma_upper = (
        harmonic_upper
        - 11 * log2_lower
        - Q // (2 * N)
        + ceil_div(Q, 12 * N * N)
    )

    # Thus gamma < 0.577215664902, with an integer-only check.
    if not gamma_upper * 10**12 < 577_215_664_902 * Q:
        raise ArithmeticError("Euler-constant bound could not be certified")

    x = Fraction(577_215_664_902, 10**12)
    t = Fraction(1)
    upper = t
    for k in range(1, 21):
        t *= x / k
        upper += t

    # Geometric upper bound for sum_{k=21}^infty x^k/k!.
    remainder = t * (x / 21) / (1 - x / 22)
    if not upper + remainder < Fraction(EGAMMA_NUM, EGAMMA_DEN):
        raise ArithmeticError("Exponential upper bound could not be certified")


def prime_sieve(limit):
    """Return a primality table and all primes from 13 to limit."""
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if is_prime[p]:
            is_prime[p * p::p] = b"\x00" * ((limit - p * p) // p + 1)
    return is_prime, [p for p in range(A + 1, limit + 1) if is_prime[p]]


def build_u_events(primes):
    """Build U_hat exactly in structure, rounded upwards term by term.

    For each p >= 13 and k >= 1, let m=(k+1)*p.
    - The generic contribution c_k/p^2 begins when n >= m+2,
      with c_1=1/2, c_k=1/k for k>=2.
    - The two exceptional sums contribute 1/(k*p^2) at n=m and n=m+1.

    Thus the exact U_hat(n) is the cumulative sum of generic events,
    plus the exceptional events at that n. We round each contribution
    UP to a multiple of 1/Q, ensuring an upper bound for U_hat(n).
    """
    generic_events = array("Q", [0]) * (N_END + 1)
    exceptional_events = array("Q", [0]) * (N_END + 1)

    for p in primes:
        if 2 * p > N_END:
            break
        p2 = p * p
        for k in range(1, N_END // p):
            m = (k + 1) * p
            exceptional = ceil_div(Q, k * p2)
            exceptional_events[m] += exceptional
            if m + 1 <= N_END:
                exceptional_events[m + 1] += exceptional
            if m + 2 <= N_END:
                generic = ceil_div(Q, (2 if k == 1 else k) * p2)
                generic_events[m + 2] += generic

    return generic_events, exceptional_events


def exp_mu_lower(mu_lower):
    """Rigorous lower bound for Q*exp(mu), where mu >= mu_lower/Q."""
    term = Q
    total = Q
    for k in range(1, 21):
        term = (term * mu_lower) // (k * Q)
        total += term
    return total


def verifies_threshold(u_upper, exp_mu_scaled):
    """Check 1 - e^gamma/e^mu - U_hat > TARGET_NUM/TARGET_DEN.

    u_upper >= Q*U_hat, exp_mu_scaled <= Q*exp(mu), and
    EGAMMA_NUM/EGAMMA_DEN > exp(gamma). All operations are integers.
    """
    left = ((TARGET_DEN - TARGET_NUM) * Q - TARGET_DEN * u_upper)
    return (
        left > 0
        and left * EGAMMA_DEN * exp_mu_scaled
        > TARGET_DEN * EGAMMA_NUM * Q * Q
    )


def main():
    check_euler_constant_bound()
    print("Certified bound: exp(gamma) < 1.781072418", flush=True)
    print("Building prime and cycle-count tables...", flush=True)
    is_prime, primes = prime_sieve(N_END)
    generic, exceptional = build_u_events(primes)

    # At n=N_START, sum all primes p <= N_START-3 in mu.
    # Terms are rounded DOWN to give a lower bound for mu.
    mu_lower = sum(Q // p for p in primes if p <= N_START - 3)
    generic_upper = sum(generic[:N_START + 1])
    exp_lower = exp_mu_lower(mu_lower)

    min_approx = float("inf")
    min_n = None

    for n in range(N_START, N_END + 1):
        if n > N_START:
            generic_upper += generic[n]
            if is_prime[n - 3] and n - 3 > A:
                mu_lower += Q // (n - 3)
                exp_lower = exp_mu_lower(mu_lower)

        u_upper = generic_upper + exceptional[n]
        if not verifies_threshold(u_upper, exp_lower):
            print("NOT CERTIFIED at n =", n)
            sys.exit(1)

        # For reporting the approximate minimum ONLY. This is not part
        # of the rigorous certificate, which uses the integer test above.
        approximate_L = 1.0 - exp(EULER_GAMMA_APPROX - mu_lower / Q) - u_upper / Q
        if approximate_L < min_approx:
            min_approx = approximate_L
            min_n = n

        if (n - N_START + 1) % PROGRESS_EVERY == 0:
            print(f"Certified through n={n:,}", flush=True)

    print()
    print("SUCCESS: every n from 9053 through 1199999 is certified.")
    print("For every n, L(n;12,n-3) > 0.333777 > 1/3.")
    print(f"Approximate minimum: {min_approx:.12f} at n={min_n}.")
    print("All certification comparisons used integer arithmetic.")


if __name__ == "__main__":
    main()
