
# prime-plus-nonsquarefree

Computational verification accompanying the paper
*Prime plus a non-square-free integer: an elementary approach*.

The paper proves that every integer $n > 24$ can be written as
the sum of a prime and a positive integer that is not square-free.

## Verification scripts

- `compute_small_prime_square_values.py`:
  Computes $P(q^2)$ and $q^-\#$ for $q = 5, 7, 11, 13$.
- `verify_small_representations.py`:
  Verifies the required representations for all integers
  $24 < n \leq 3617$ satisfying $6 \mid n$.
- `verify_prime_square_bound.py`:
  Verifies $P(q^2) < q^-\#$ for every prime
  $17 \leq q \leq 14\,347$.

## Running the scripts

The first two scripts require only the Python standard library.
The third requires NumPy and Numba, which can be installed using:

    python -m pip install numpy numba

Each script can then be run using:

    python compute_small_prime_square_values.py
    python verify_small_representations.py
    python verify_prime_square_bound.py

The third script creates a new CSV file and will not overwrite
an existing results file.

## Runtime benchmark

`verify_prime_square_bound.py` completed the verification of all
1,677 primes in 5.634 hours on:

- **CPU:** AMD Ryzen 5 5600X (3.70 GHz base clock)
- **Execution:** single-threaded
- **OS:** Windows 11, 64-bit

The other two scripts each completed in under one second.

## Computational results

The complete results of `verify_prime_square_bound.py` are
included in `prime_square_bound_results.csv`.

All 1,677 primorial inequalities were verified successfully.
The largest computed value of $P(q^2)$ was 94,121,869,957,
attained at $q = 14\,293$.
