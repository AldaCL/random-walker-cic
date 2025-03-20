// This code is an adaptation by @aldaCortes to work with Apple silicon Clang due to the
// lack of _unit64_t in the original code on the paper.
// Original code paper https://www.hvks.com/Numerical/Downloads/HVE%20Practical%20implementation%20of%20Spigot%20Algorithms%20for%20transcendental%20constants.pdf 

// This code uses spigot algorithm to compute the first n digits of PI sequentially. 
// The algorithm is based on the paper "Practical implementation of Spigot Algorithms for transcendental constants" by H. V. Koning and J. W. van Holten.



#include <string>
#include <iostream>
#include <climits>
#include <cstdio>
#include <cstdint>

std::string pi_spigot_64(const int digits, int no_dig = 4)
{
    static uint64_t f_table[9] = {0, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    static uint64_t f2_table[9] = {0, 2, 20, 200, 2000, 20000, 200000, 2000000, 20000000};
    const int TERMS = (10 * no_dig / 3 + 1);
    bool first_time = true;
    bool overflow_flag = false;
    char buffer[32];
    std::string ss;
    long b, c;
    int carry, no_carry = 0;
    uint64_t f, f2;
    uint64_t dig_n = 0;
    uint64_t e = 0;
    uint64_t acc = 0;
    uint64_t g = 0;
    uint64_t tmp64;

    ss.reserve(digits + 16);
    if (no_dig > 8) no_dig = 8;
    if (no_dig < 1) no_dig = 1;
    
    c = (digits / no_dig + 1) * no_dig;
    if (no_dig == 1) c++;
    c = (c / no_dig + 1) * TERMS;

    f = f_table[no_dig];
    f2 = f2_table[no_dig];

    uint64_t *a = new uint64_t[c];

    for (; (b = c -= TERMS) > 0 && !overflow_flag; first_time = false)
    {
        for (; --b > 0 && !overflow_flag;)
        {
            if (acc > ULLONG_MAX / b) overflow_flag = true;
            acc *= b;
            tmp64 = f;
            if (first_time)
                tmp64 *= f2;
            else
                tmp64 *= a[b];

            if (acc > ULLONG_MAX - tmp64) overflow_flag = true;
            acc += tmp64;
            g = b + b - 1;
            a[b] = acc % g;
            acc /= g;
        }

        dig_n = (uint64_t)(e + acc / f);
        carry = (unsigned)(dig_n / f);
        dig_n %= f;

        if (carry > 0)
        {
            ++no_carry;
            for (size_t i = ss.length(); carry > 0 && i > 0; --i)
            {
                int new_digit = (ss[i - 1] - '0') + carry;
                carry = new_digit / 10;
                ss[i - 1] = new_digit % 10 + '0';
            }
        }

        snprintf(buffer, sizeof(buffer), "%0*llu", no_dig, static_cast<unsigned long long>(dig_n));        ss += std::string(buffer);
        if (first_time) ss.insert(1, ".");

        acc = acc % f;
        e = (uint64_t)acc;
    }

    ss.erase(digits + 1);
    if (overflow_flag) ss = std::string("Overflow:") + ss;

    delete[] a;
    return ss;
}


int main()
{
    int digits = 50; // Example: Compute first 50 digits of PI
    std::string pi = pi_spigot_64(digits);
    std::cout << "PI (first " << digits << " digits): " << pi << std::endl;
    return 0;
}