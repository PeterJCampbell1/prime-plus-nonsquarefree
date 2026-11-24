"""Compute P(q^2) and compare it with the product of the primes below q.

For each prime q, use a segmented sieve to find prime representatives
of the reduced residue classes modulo q^2. Once at most 16,384 classes
remain, search each remaining class directly using exact trial division
to find its least prime representative.

Run this file after installing numpy and numba. A new CSV is created for
each run; existing files are never overwritten or automatically reused.
"""

import csv
import math
import time

import numpy as np
from numba import njit


MIN_Q = 17
MAX_Q = 14_347
OUTPUT_FILE = "prime_square_bound_results1.csv"

SEARCH_CAP = 10**12         # Largest prime candidate allowed, inclusive.
SEGMENT_ODDS = 524_288      # Number of odd integers stored in each sieve segment.
CHUNK_SPAN = 1_073_741_824  # Maximum interval before checking to switch to trial division.
TAIL_CLASSES = 16_384       # Switch to trial division at this many missing classes.


def primes_less_than(limit):
    """Return a list of all primes p < limit using a basic sieve."""
    if limit <= 2:
        return []

    is_prime = bytearray([1]) * limit
    is_prime[0] = is_prime[1] = 0

    for d in range(2, math.isqrt(limit - 1) + 1):
        if is_prime[d]:
            count = (limit - 1 - d * d) // d + 1
            is_prime[d * d : limit : d] = bytearray(count)

    return [n for n in range(2, limit) if is_prime[n]]  


@njit(cache=True)
def sieve_interval(q, covered, missing, start, stop, small_primes, segment_odds):
    """Sieve the odd integers from start through stop, in increasing order.

    start is odd and at least 3. covered[a] is 1 when a is not reduced or
    already has a prime representative. Return the next odd integer,
    the number of missing classes, and P(q^2) if coverage is complete.
    """
    modulus = q * q

    while start <= stop:
        length = min(segment_odds, (stop - start) // 2 + 1)
        end = start + 2 * (length - 1)
        composite = np.zeros(length, dtype=np.uint8)

        # Entry i represents the odd integer start + 2*i.
        for p in small_primes:
            if p == 2:
                continue
            if p * p > end:
                break

            first = max(p * p, ((start + p - 1) // p) * p)
            if first % 2 == 0:
                first += p

            for i in range((first - start) // 2, length, p):
                composite[i] = 1

        for i in range(length):
            if composite[i] == 0:
                p = start + 2 * i
                residue = p % modulus

                if covered[residue] == 0:
                    covered[residue] = 1
                    missing -= 1
                    if missing == 0:
                        return p + 2, 0, p

        start = end + 2

    return start, missing, -1


@njit(cache=True)
def is_prime_by_trial_division(n, small_primes):
    """Test n exactly; small_primes must include all primes up to sqrt(n)."""
    if n < 2:
        return False
    for p in small_primes:
        if p * p > n:
            break
        if n % p == 0:
            return False
    return True


@njit(cache=True)
def finish_remaining_classes(q, covered, missing, start, search_cap, small_primes):
    """Search each uncovered progression, starting beyond the sieved range.

    All primes below start have already been processed. Thus the first prime
    found in each remaining class is its least representative, and exceeds
    all representatives found in the sieve stage. Return -1 on an incomplete
    search; missing then counts classes still without a verified witness.
    """
    modulus = q * q
    largest_prime = -1

    for residue in range(modulus):
        if covered[residue]:
            continue

        # The first integer >= start congruent to residue modulo q^2.
        candidate = start + (residue - start) % modulus
        step = modulus
        if q != 2:
            if candidate % 2 == 0:
                candidate += modulus
            step = 2 * modulus   # Skip even candidates; 2 was handled earlier.

        while candidate <= search_cap:
            if is_prime_by_trial_division(candidate, small_primes):
                covered[residue] = 1
                missing -= 1
                largest_prime = max(largest_prime, candidate)
                break
            candidate += step
        else:
            return -1, missing

    return largest_prime, missing



def largest_least_prime(q, small_primes):
    """Return (P(q^2), 0), or (-1, missing) if SEARCH_CAP is reached.

    q must be prime. small_primes must contain every prime <= sqrt(SEARCH_CAP).
    """
    modulus = q * q
    covered = np.zeros(modulus, dtype=np.uint8)
    covered[::q] = 1       # Exclude non-reduced residue classes.
    missing = q * (q - 1)

    if q != 2:            # Account for the only even prime separately.
        covered[2] = 1
        missing -= 1

    start = 3

    while start <= SEARCH_CAP:
        # Check for the sparse tail regularly, even for small q.
        stop = min(SEARCH_CAP, start + min(CHUNK_SPAN, 8 * modulus) - 1)

        start, missing, result = sieve_interval(
            q, covered, missing, start, stop, small_primes, SEGMENT_ODDS
        )

        if missing == 0:
            return int(result), 0

        if 0 < missing <= TAIL_CLASSES and start <= SEARCH_CAP:
            result, missing = finish_remaining_classes(
                q, covered, missing, start, SEARCH_CAP, small_primes
            )
            return int(result), int(missing)

    return -1, int(missing)


def warm_up(small_primes):
    """Compile both numerical stages before timing, and check P(25) = 149."""
    covered = np.zeros(25, dtype=np.uint8)
    covered[::5] = 1
    covered[2] = 1
    start, missing, _ = sieve_interval(5, covered, 19, 3, 31, small_primes, 16)
    result, missing = finish_remaining_classes(5, covered, missing, start, 1000, small_primes)
    if (result, missing) != (149, 0):
        raise RuntimeError("The P(25) self-check failed.")


def main():
    if not 2 <= MIN_Q <= MAX_Q:
        raise ValueError("Require 2 <= MIN_Q <= MAX_Q.")
    # Keep every integer in the compiled loops safely within signed int64.
    if not 2 <= SEARCH_CAP <= 2**40:
        raise ValueError("Keep SEARCH_CAP between 2 and 2**40.")
    if MAX_Q * MAX_Q > 1024**3:
        raise ValueError("The coverage array would exceed 1 GiB; review the range first.")
    if not 1 <= SEGMENT_ODDS <= 2**26 or not 2 <= CHUNK_SPAN <= 2**40 or TAIL_CLASSES < 0:
        raise ValueError("Invalid sieve or tail settings.")

    q_values = primes_less_than(MAX_Q + 1)
    if q_values[-1] < MIN_Q:
        raise ValueError("There are no primes in the requested range.")

    # Enough primes for every candidate, and for the small warm-up check.
    small_primes = np.array(
        primes_less_than(math.isqrt(max(SEARCH_CAP, 1000)) + 1), dtype=np.int64
    )

    previous_primorial = 1
    max_P = -1
    max_q = None
    total_seconds = 0
    checked = 0
    all_hold = True

    # 'x' creates a fresh file and refuses to overwrite an existing result.
    with open(OUTPUT_FILE, "x", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow([
            "q", "P_q_squared", "extremal_residue", "P_lt_primorial",
            "seconds", "status", "missing", "search_cap"
        ])
        output.flush()  # Write the CSV header to the file immediately.

        # Compile the Numba functions and check them on a small example
        # before timing the main computations.
        print("Preparing the compiled loops...", flush=True)
        warm_up(small_primes)

        for q in q_values:
            q_minus_primorial = previous_primorial
            previous_primorial *= q
            if q < MIN_Q:
                continue

            print(f"Starting q={q}; {q * (q - 1):,} reduced classes", flush=True)
            started = time.perf_counter()
            P_q, missing = largest_least_prime(q, small_primes)
            elapsed = time.perf_counter() - started
            total_seconds += elapsed

            if missing:
                writer.writerow([q, "", "", "", f"{elapsed:.6f}", "INCOMPLETE", missing, SEARCH_CAP])
                output.flush()
                raise SystemExit(
                    f"Search cap reached for q={q}. This q is INCOMPLETE, not a counterexample.\n"
                    "Completed rows are saved. Increase SEARCH_CAP before retrying."
                )

            extremal_class = P_q % (q * q)
            inequality_holds = P_q < q_minus_primorial
            writer.writerow([
                q, P_q, extremal_class, inequality_holds, f"{elapsed:.6f}",
                "COMPLETE", 0, SEARCH_CAP
            ])
            output.flush()  # Save each completed q, rather than waiting for the entire run.
            print(
                f"q={q}: P(q^2)={P_q:,}; extremal class={extremal_class}; "
                f"P(q^2)<q^-# is {inequality_holds}; elapsed {elapsed:.3f}s",
                flush=True
            )

            checked += 1
            if not inequality_holds:
                all_hold = False
            if P_q > max_P:
                max_P = P_q
                max_q = q

    print(f"\nCompleted {checked} prime values in {total_seconds / 3600:.3f} hours of computation.")
    print(f"Maximum P(q^2) = {max_P:,}, attained at q = {max_q}.")
    print(f"All primorial inequalities hold: {all_hold}")
    print(f"Results saved to {OUTPUT_FILE}")
    if not all_hold:
        raise SystemExit("At least one least-prime inequality failed; inspect the CSV.")


if __name__ == "__main__":
    try:
        main()
    except FileExistsError:
        raise SystemExit(f"{OUTPUT_FILE} already exists. Rename it or change OUTPUT_FILE; nothing was overwritten.")
    except KeyboardInterrupt:
        raise SystemExit("Interrupted. Completed CSV rows are saved; the unfinished q has no result.")
