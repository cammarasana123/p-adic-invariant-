# A p-adic invariant of primes via the Möbius transform of the Lucas sequence**


This repository contains the Python script used for the computational verification in the paper

**A p-adic invariant of primes via the Möbius transform of the Lucas sequence**  

The script verifies, over a finite numerical range, the computational statements reported in the paper:

1. **Main Theorem, part (1):** convergence and stabilization modulo \(p^n\);
2. **Main Theorem, part (2):** exact rate of convergence;
3. **Digit Law:** closed formula for the incremental digits and the split/inert pattern.

The computations use modular fast doubling for Fibonacci and Lucas numbers. The script uses only the Python standard library; no external packages are required.

## Requirements

Python 3.9 or newer is recommended.

No additional libraries are needed.

You can check your Python version with:

```bash
python3 --version
```

## Files

```text
cammarasana_reproducibility.py
    Python script for the computational verification.

cammarasana_reproducibility_output.txt
    Output of the default run.
```

## Default verification

Run the default verification with:

```bash
python3 cammarasana_reproducibility.py
```

The default range is:

```text
7 <= p <= 100000
1 <= n <= 7
1 <= k <= 7
```

The default run tests all primes \(p\) in the range \(7 \le p \le 100000\), excluding \(p=2,3,5\).

The output begins with the following summary:

```text
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
```

A sample table is printed after the summary unless it is suppressed with `--no-table`.

## Command-line options

Display the available options with:

```bash
python3 cammarasana_reproducibility.py --help
```

Available options:

```text
--max-prime N
    Largest prime tested.
    Default: 100000

--max-k K
    Largest k tested.
    Default: 7

--max-n N
    Largest n used for stabilization modulo p^n.
    Default: 7

--samples LIST
    Comma-separated list of sample primes shown in the output table.
    Default: 7,11,13,17,19,23,29,37

--sample-k K
    Number of k-values shown in the sample table.
    Default: 5

--workers W
    Number of worker processes.
    Default: 1

--no-table
    Suppress the sample table.
```

The parameter `--max-n` must not exceed `--max-k`.

## Examples

Run the default verification:

```bash
python3 cammarasana_reproducibility.py
```

Run the default verification without the sample table:

```bash
python3 cammarasana_reproducibility.py --no-table
```

Run a smaller test:

```bash
python3 cammarasana_reproducibility.py --max-prime 1000 --max-k 5 --max-n 5
```

Run the default range using four worker processes:

```bash
python3 cammarasana_reproducibility.py --workers 4
```

Change the sample primes shown in the output table:

```bash
python3 cammarasana_reproducibility.py --samples 7,11,13,37,41
```

## Reproducibility statement

The default run verifies the computational claims reported in the paper for all primes \(7 \le p \le 100000\), with \(1 \le n \le 7\) and \(1 \le k \le 7\).

The reported default run gives zero mismatches in all tests.
