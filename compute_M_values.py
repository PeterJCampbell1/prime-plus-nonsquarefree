MIN_Q = 2
MAX_Q = 500


def primes_less_than(limit):
    """Return a list of all primes p < limit using a basic sieve."""
    if limit <= 2:
        return []

    is_prime = bytearray([1]) * limit
    is_prime[0] = is_prime[1] = 0

    for d in range(2, int(limit ** 0.5) + 1):
        if is_prime[d]:
            is_prime[d * d : limit : d] = bytearray(
                len(is_prime[d * d : limit : d])
            )

    return [n for n in range(2, limit) if is_prime[n]]


def least_prime_representatives(q):
    """Return the least prime in every reduced residue class modulo q^2."""
    modulus = q * q
    number_of_classes = q * (q - 1)
    search_limit = max(100, 2 * modulus)

    while True:
        representatives = {}

        for p in primes_less_than(search_limit):
            residue = p % modulus

            if residue % q != 0 and residue not in representatives:
                representatives[residue] = p

                if len(representatives) == number_of_classes:
                    return representatives

        search_limit *= 2


def main():
    q_values = primes_less_than(MAX_Q + 1)
    previous_primorial = 1

    print("q,M(q),q^-#,residue class attaining M(q),M(q) < q^-#")

    for q in q_values:
        q_minus_primorial = previous_primorial
        previous_primorial *= q

        if q < MIN_Q:
            continue

        representatives = least_prime_representatives(q)
        M_q = max(representatives.values())
        extremal_class = M_q % (q * q)

        inequality_holds = M_q < q_minus_primorial

        print(
            f"{q},{M_q},{q_minus_primorial},{extremal_class},"
            f"{inequality_holds}"
        )


if __name__ == "__main__":
    main()
