LOWER_BOUND = 24
UPPER_BOUND = 3617


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


def least_square_prime_divisor(number, primes):
    """Return the least prime ell for which ell^2 divides number, if one exists."""
    for ell in primes:
        square = ell * ell

        if square > number:
            return None

        if number % square == 0:
            return ell

    return None


def find_solution(n, primes):
    """Find (p, r, ell) with n = p + r and ell^2 dividing r."""
    for p in primes:
        if p >= n:
            break

        r = n - p
        ell = least_square_prime_divisor(r, primes)

        if ell is not None:
            return p, r, ell

    return None


def main():
    primes = primes_less_than(UPPER_BOUND)

    numbers_checked = 0

    for n in range(LOWER_BOUND + 1, UPPER_BOUND + 1):
        if n % 6 != 0:
            continue

        solution = find_solution(n, primes)

        if solution is None:
            raise RuntimeError(
                f"No representation was found for n = {n}."
            )

        p, r, ell = solution

        print(f"{n}={p}+{r}, where {r} is divisible by {ell}^2 = {ell**2}.")
        numbers_checked += 1

    print()
    print(
        f"Verified all {number_checked} integers n satisfying "
        f"{LOWER_BOUND} < n <= {UPPER_BOUND} and 6 | n."
    )


if __name__ == "__main__":
    main()