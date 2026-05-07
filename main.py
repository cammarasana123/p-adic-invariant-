#!/usr/bin/env python3
"""
Computational verification for
"A p-adic invariant of primes via the Mobius transform of the Lucas sequence".

For every prime p in the chosen range, p >= 7, the script verifies:

  Main Theorem, part (1): p^n-stabilisation of w_k(p), equivalently
      u_k(p) == chi(p)^(k-1) * lambda(p) mod p^n
    in the tested range 1 <= n <= N and n <= k <= K.

  Main Theorem, part (2): the exact convergence rate
      v_p(w_{k+1}(p) - w_k(p)) = k
    in the tested range 1 <= k <= K.

  Digit Law:
      d_k(p) == chi(p)^(k-1) * c_2(p) mod p,
      c_2(p) = d_1(p) == (5/8) * (F_{p-chi(p)} / p)^2 mod p,
    in the tested range 1 <= k <= K.

Default range:
    all primes 7 <= p <= 100000,
    Main Theorem exact rate and Digit Law: 1 <= k <= 7,
    Main Theorem congruence form: 1 <= n <= 7 and n <= k <= 7.

The computations use modular fast doubling for Fibonacci and Lucas numbers.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

def fib_pair_mod(n: int, m: int) -> Tuple[int, int]:
    """Return (F_n, F_{n+1}) modulo m by fast doubling."""
    if n == 0:
        return 0, 1 % m

    a, b = fib_pair_mod(n >> 1, m)
    c = a * (2 * b - a) % m
    d = (a * a + b * b) % m

    if n & 1:
        return d, (c + d) % m
    return c, d


def fib_mod(n: int, m: int) -> int:
    """Return F_n modulo m."""
    return fib_pair_mod(n, m)[0]


def lucas_mod(n: int, m: int) -> int:
    """Return L_n modulo m, using L_n = 2 F_{n+1} - F_n."""
    fn, fn1 = fib_pair_mod(n, m)
    return (2 * fn1 - fn) % m


def legendre_5(p: int) -> int:
    """Return the Legendre symbol (5/p), for odd primes p != 5."""
    r = pow(5, (p - 1) // 2, p)
    return -1 if r == p - 1 else r


def chi_power(chi: int, exponent: int) -> int:
    """Return chi^exponent as the integer +1 or -1."""
    if chi == 1 or exponent % 2 == 0:
        return 1
    return -1


def primes_up_to(n: int) -> List[int]:
    """Return all primes <= n."""
    if n < 2:
        return []

    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(n**0.5)
    for q in range(2, limit + 1):
        if sieve[q]:
            start = q * q
            sieve[start : n + 1 : q] = b"\x00" * (((n - start) // q) + 1)

    return [q for q in range(2, n + 1) if sieve[q]]


def vp_residue(x: int, p: int, modulus_exponent: int) -> int:
    """
    Return v_p(x) for a residue modulo p^modulus_exponent.
    If x is 0 modulo p^modulus_exponent, return modulus_exponent.
    """
    if x == 0:
        return modulus_exponent

    v = 0
    while v < modulus_exponent and x % p == 0:
        x //= p
        v += 1
    return v


# ---------------------------------------------------------------------------
# Definitions from the paper
# ---------------------------------------------------------------------------

def normalized_u_mod(p: int, k: int, precision: int) -> int:
    """Return u_k(p) = B(p^k) / p^k modulo p^precision."""
    modulus = p ** (k + precision)
    pk = p**k

    b_mod = (lucas_mod(pk, modulus) - lucas_mod(pk // p, modulus)) % modulus

    if b_mod % pk != 0:
        raise ArithmeticError(f"B(p^k) is not divisible by p^k for p={p}, k={k}")

    return (b_mod // pk) % (p**precision)


def twisted_w_mod(p: int, k: int, chi: int, precision: int) -> int:
    """Return w_k(p) = u_k(p) * chi(p)^(1-k) modulo p^precision."""
    u = normalized_u_mod(p, k, precision)
    return (u * chi_power(chi, 1 - k)) % (p**precision)


def c2_formula(p: int, chi: int) -> int:
    """
    Return c_2(p) = (5/8) * (F_{p-chi(p)} / p)^2 modulo p.
    """
    idx = p - chi
    f = fib_mod(idx, p * p)

    if f % p != 0:
        raise ArithmeticError(f"F_{{p-chi}} is not divisible by p for p={p}")

    quotient = (f // p) % p
    inv8 = pow(8, p - 2, p)
    return (5 * inv8 * quotient * quotient) % p


def digit_formula(p: int, k: int, chi: int, c2: int) -> int:
    """Return chi(p)^(k-1) * c_2(p) modulo p."""
    return (chi_power(chi, k - 1) * c2) % p


class WCache:
    """Cache values w_k(p) modulo a fixed power p^precision for a fixed p."""

    def __init__(self, p: int, chi: int, precision: int) -> None:
        self.p = p
        self.chi = chi
        self.precision = precision
        self.modulus = p**precision
        self.values: Dict[int, int] = {}

    def w(self, k: int) -> int:
        if k not in self.values:
            self.values[k] = twisted_w_mod(self.p, k, self.chi, self.precision)
        return self.values[k]


def digit_and_valuation(p: int, k: int, cache: WCache) -> Tuple[Optional[int], int]:
    """Return (d_k(p), v_p(w_{k+1}-w_k))."""
    diff = (cache.w(k + 1) - cache.w(k)) % cache.modulus
    valuation = vp_residue(diff, p, cache.precision)

    pk = p**k
    if diff % pk != 0:
        return None, valuation

    return (diff // pk) % p, valuation


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

@dataclass
class SampleRow:
    p: int
    chi: int
    lambda_mod_p: int
    c2: int
    valuations: List[int]
    digits: List[int]


@dataclass
class VerificationSummary:
    max_prime: int
    max_k: int
    max_n: int
    primes_tested: int
    split_primes: int
    inert_primes: int
    main_part1_cases: int
    main_part1_failures: int
    main_part2_cases: int
    main_part2_failures: int
    digit_law_cases: int
    digit_law_mismatches: int
    split_pattern_failures: int
    inert_pattern_failures: int
    c2_formula_mismatches: int


@dataclass
class PrimeResult:
    p: int
    chi: int
    main_part1_failures: int
    main_part2_failures: int
    digit_law_mismatches: int
    split_pattern_failures: int
    inert_pattern_failures: int
    c2_formula_mismatches: int
    sample_row: Optional[SampleRow]


def verify_one_prime(args: Tuple[int, int, int, Tuple[int, ...], int]) -> PrimeResult:
    p, max_k, max_n, sample_tuple, sample_k = args
    sample_set = set(sample_tuple)

    chi = legendre_5(p)
    if chi not in (-1, 1):
        raise ArithmeticError(f"Unexpected chi={chi} for p={p}")

    cache = WCache(p, chi, precision=max_k + 1)
    c2 = c2_formula(p, chi)

    main_part1_failures = 0
    main_part2_failures = 0
    digit_law_mismatches = 0
    split_pattern_failures = 0
    inert_pattern_failures = 0
    c2_formula_mismatches = 0

    for n in range(1, max_n + 1):
        pn = p**n
        lambda_mod_pn = cache.w(n) % pn
        for k in range(n, max_k + 1):
            if cache.w(k) % pn != lambda_mod_pn:
                main_part1_failures += 1

    data_digits: List[int] = []
    valuations: List[int] = []

    for k in range(1, max_k + 1):
        digit, valuation = digit_and_valuation(p, k, cache)
        expected_digit = digit_formula(p, k, chi, c2)

        valuations.append(valuation)
        data_digits.append(digit if digit is not None else -1)

        if valuation != k:
            main_part2_failures += 1
        if digit != expected_digit:
            digit_law_mismatches += 1

    if data_digits and data_digits[0] != c2:
        c2_formula_mismatches += 1

    if chi == 1:
        for d in data_digits:
            if d != c2:
                split_pattern_failures += 1
    else:
        for k, d in enumerate(data_digits, start=1):
            expected = c2 if k % 2 == 1 else (-c2) % p
            if d != expected:
                inert_pattern_failures += 1

    sample_row = None
    if p in sample_set:
        shown_k = min(sample_k, max_k)
        sample_row = SampleRow(
            p=p,
            chi=chi,
            lambda_mod_p=cache.w(1) % p,
            c2=c2,
            valuations=valuations[:shown_k],
            digits=data_digits[:shown_k],
        )

    return PrimeResult(
        p=p,
        chi=chi,
        main_part1_failures=main_part1_failures,
        main_part2_failures=main_part2_failures,
        digit_law_mismatches=digit_law_mismatches,
        split_pattern_failures=split_pattern_failures,
        inert_pattern_failures=inert_pattern_failures,
        c2_formula_mismatches=c2_formula_mismatches,
        sample_row=sample_row,
    )


def verify(
    max_prime: int,
    max_k: int,
    max_n: int,
    sample_primes: Iterable[int],
    sample_k: int,
    workers: int = 1,
) -> Tuple[VerificationSummary, List[SampleRow]]:
    if max_prime < 7:
        raise ValueError("max_prime must be at least 7")
    if max_k < 1:
        raise ValueError("max_k must be at least 1")
    if max_n < 1:
        raise ValueError("max_n must be at least 1")
    if max_n > max_k:
        raise ValueError("max_n must not exceed max_k")

    primes = [p for p in primes_up_to(max_prime) if p not in (2, 3, 5)]
    sample_tuple = tuple(sample_primes)

    main_part1_cases_per_prime = sum(max_k - n + 1 for n in range(1, max_n + 1))
    jobs = [(p, max_k, max_n, sample_tuple, sample_k) for p in primes]

    if workers and workers > 1:
        with mp.Pool(processes=workers) as pool:
            results = list(pool.imap_unordered(verify_one_prime, jobs, chunksize=32))
    else:
        results = [verify_one_prime(job) for job in jobs]

    split_primes = sum(1 for r in results if r.chi == 1)
    inert_primes = sum(1 for r in results if r.chi == -1)
    main_part1_failures_count = sum(r.main_part1_failures for r in results)
    main_part2_failures_count = sum(r.main_part2_failures for r in results)
    digit_law_mismatches_count = sum(r.digit_law_mismatches for r in results)
    split_pattern_failures_count = sum(r.split_pattern_failures for r in results)
    inert_pattern_failures_count = sum(r.inert_pattern_failures for r in results)
    c2_formula_mismatches_count = sum(r.c2_formula_mismatches for r in results)
    sample_rows = [r.sample_row for r in results if r.sample_row is not None]

    summary = VerificationSummary(
        max_prime=max_prime,
        max_k=max_k,
        max_n=max_n,
        primes_tested=len(primes),
        split_primes=split_primes,
        inert_primes=inert_primes,
        main_part1_cases=len(primes) * main_part1_cases_per_prime,
        main_part1_failures=main_part1_failures_count,
        main_part2_cases=len(primes) * max_k,
        main_part2_failures=main_part2_failures_count,
        digit_law_cases=len(primes) * max_k,
        digit_law_mismatches=digit_law_mismatches_count,
        split_pattern_failures=split_pattern_failures_count,
        inert_pattern_failures=inert_pattern_failures_count,
        c2_formula_mismatches=c2_formula_mismatches_count,
    )

    return summary, sorted(sample_rows, key=lambda row: row.p)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def print_sample_table(rows: List[SampleRow]) -> None:
    if not rows:
        return

    shown_k = len(rows[0].digits)
    print("Sample table")
    print("-" * 96)
    header = (
        f"{'p':>7} {'chi':>5} {'lambda mod p':>13} {'c_2':>7} | "
        + " ".join(f"v_{k}".rjust(5) for k in range(1, shown_k + 1))
        + " | "
        + " ".join(f"d_{k}".rjust(5) for k in range(1, shown_k + 1))
    )
    print(header)
    print("-" * len(header))

    for row in rows:
        valuations = " ".join(f"{v:>5}" for v in row.valuations)
        digits = " ".join(f"{d:>5}" for d in row.digits)
        print(
            f"{row.p:>7} {row.chi:>5} {row.lambda_mod_p:>13} {row.c2:>7} | "
            f"{valuations} | {digits}"
        )

    print()


def print_summary(summary: VerificationSummary) -> None:
    part1_verified = summary.main_part1_failures == 0
    part2_verified = summary.main_part2_failures == 0
    digit_law_verified = (
        summary.digit_law_mismatches == 0
        and summary.c2_formula_mismatches == 0
        and summary.split_pattern_failures == 0
        and summary.inert_pattern_failures == 0
    )
    passed = part1_verified and part2_verified and digit_law_verified

    print("Computational verification")
    print("=" * 96)
    print(
        f"Range: 7 <= p <= {summary.max_prime}, "
        f"1 <= n <= {summary.max_n}, "
        f"1 <= k <= {summary.max_k}"
    )
    print(f"Primes tested: {summary.primes_tested}")
    print()

    print(
        "Main Theorem, part (1), convergence:               "
        f"{'VERIFIED' if part1_verified else 'FAILED'}"
    )
    print(f"  Stabilization mismatches modulo p^n:             {summary.main_part1_failures}")
    print()

    print(
        "Main Theorem, part (2), exact rate of convergence: "
        f"{'VERIFIED' if part2_verified else 'FAILED'}"
    )
    print(f"  Valuation mismatches:                            {summary.main_part2_failures}")
    print()

    print(f"Digit Law:                                        {'VERIFIED' if digit_law_verified else 'FAILED'}")
    print(f"  Digit formula mismatches:                        {summary.digit_law_mismatches}")
    print(f"  Closed formula mismatches for c_2:               {summary.c2_formula_mismatches}")
    print(f"  Split-pattern failures:                          {summary.split_pattern_failures}")
    print(f"  Inert-pattern failures:                          {summary.inert_pattern_failures}")
    print()

    print(f"Final status:                                     {'PASSED' if passed else 'FAILED'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify the Main Theorem and Digit Law for the Lucas Mobius-transform paper."
    )
    parser.add_argument("--max-prime", type=int, default=100000, help="largest prime tested; default: 100000")
    parser.add_argument("--max-k", type=int, default=7, help="largest k tested; default: 7")
    parser.add_argument("--max-n", type=int, default=7, help="largest n for p^n-stabilisation; default: 7")
    parser.add_argument(
        "--samples",
        type=str,
        default="7,11,13,17,19,23,29,37",
        help="comma-separated sample primes shown in the output table",
    )
    parser.add_argument("--sample-k", type=int, default=5, help="number of k-values shown in the sample table")
    parser.add_argument("--workers", type=int, default=1, help="number of worker processes; default: 1")
    parser.add_argument("--no-table", action="store_true", help="suppress the sample table")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_primes = [int(x) for x in args.samples.split(",") if x.strip()]

    summary, sample_rows = verify(
        max_prime=args.max_prime,
        max_k=args.max_k,
        max_n=args.max_n,
        sample_primes=sample_primes,
        sample_k=args.sample_k,
        workers=args.workers,
    )

    print_summary(summary)
    if not args.no_table:
        print()
        print_sample_table(sample_rows)


if __name__ == "__main__":
    main()
