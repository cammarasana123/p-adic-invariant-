**A p-adic invariant of primes via the Möbius transform of the Lucas sequence**

This repository contains the Python script used for the computational verification in the paper *A p-adic invariant of primes via the Möbius transform of the Lucas sequence*.

The script verifies, over a finite numerical range, the computational statements reported in the paper: Main Theorem, part (1), namely convergence and stabilization modulo p^n; Main Theorem, part (2), namely the exact rate of convergence; and the Digit Law, including the closed formula and the split/inert digit pattern.

The computations use modular fast doubling for Fibonacci and Lucas numbers. The script uses only the Python standard library; no external packages are required. Python 3.9 or newer is recommended.

To run the default verification, use: python3 cammarasana_reproducibility.py

The default range is: 7 <= p <= 100000, 1 <= n <= 7, 1 <= k <= 7.

The default run tests all primes p in the range 7 <= p <= 100000, excluding p = 2, 3, 5. It produces an output of the form:

Computational verification
================================================================================================
Range: 7 <= p <= 100000, 1 <= n <= 7, 1 <= k <= 7
Primes tested: 9589

Main Theorem, part (1), convergence:               VERIFIED
  Stabilization mismatches modulo p^n:             0

Main Theorem, part (2), exact rate of convergence: VERIFIED
  Valuation mismatches:                            0

Digit Law:                                        VERIFIED
  Digit formula mismatches:                        0
  Closed formula mismatches for c_2:               0
  Split-pattern failures:                         0
  Inert-pattern failures:                         0

Final status:                                     PASSED

A sample table is printed after the summary unless suppressed.

Command-line options are available through: python3 cammarasana_reproducibility.py --help

Available options:

--max-prime N
    Largest prime tested. Default: 100000.

--max-k K
    Largest k tested. Default: 7.

--max-n N
    Largest n used for stabilization modulo p^n. Default: 7.

--samples LIST
    Comma-separated list of sample primes shown in the output table.
    Default: 7,11,13,17,19,23,29,37.

--sample-k K
    Number of k-values shown in the sample table. Default: 5.

--workers W
    Number of worker processes. Default: 1.

--no-table
    Suppress the sample table.

The parameter --max-n must not exceed --max-k.

Examples:

python3 cammarasana_reproducibility.py
python3 cammarasana_reproducibility.py --no-table
python3 cammarasana_reproducibility.py --max-prime 1000 --max-k 5 --max-n 5
python3 cammarasana_reproducibility.py --workers 4
python3 cammarasana_reproducibility.py --samples 7,11,13,37,41
