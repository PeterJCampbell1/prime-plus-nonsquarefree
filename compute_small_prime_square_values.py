"""Compute P(q^2) and q^-# for the small-prime table in the paper.

The default range reproduces the rows for q = 5, 7, 11, and 13.
"""

import math


MIN_Q = 5
MAX_Q = 13

SEARCH_LIMIT = 5000  # Upper bound for finding prime representatives.


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



def largest_least_prime(q):
    """Return P(q^2) and a residue class attaining it; q must be prime."""
    modulus = q * q
    number_of_classes = q * (q - 1)
    covered = bytearray(modulus)
    missing = number_of_classes

    for p in primes_less_than(SEARCH_LIMIT):
        residue = p % modulus

        if residue % q != 0 and not covered[residue]:
            covered[residue] = 1
            missing -= 1

            if missing == 0:
                return p, residue

    # The search limit was reached before all reduced residue classes were covered.
    raise RuntimeError(
        f"Not all residue classes were covered for q={q}. "
        f"{missing} classes remain."
    )



def main():
    q_values = primes_less_than(MAX_Q + 1)
    previous_primorial = 1

    print("q,P(q^2),q^-#,residue class attaining P(q^2),P(q^2) < q^-#")

    for q in q_values:
        q_minus_primorial = previous_primorial
        previous_primorial *= q

        if q < MIN_Q:
            continue

        P_q, extremal_class = largest_least_prime(q)
        inequality_holds = P_q < q_minus_primorial

        print(
            f"{q},{P_q},{q_minus_primorial},{extremal_class},"
            f"{inequality_holds}"
        )


if __name__ == "__main__":
    main()

