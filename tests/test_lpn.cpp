#include <pvac/pvac.hpp>

#include <vector>
#include <random>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace pvac;

static uint8_t bitvec_dot2(const BitVec& a, const BitVec& b) {
    size_t words = (a.nbits + 63) / 64;
    uint64_t parity = 0;
    for (size_t i = 0; i < words; ++i) {
        parity ^= (uint64_t)__builtin_popcountll(a.w[i] & b.w[i]);
    }
    return (uint8_t)(parity & 1ull);
}

static BitVec random_bitvec(int m, std::mt19937_64& rng) {
    BitVec v = BitVec::make(m);
    size_t words = (m + 63) / 64;
    for (size_t i = 0; i < words; ++i) {
        v.w[i] = rng();
    }
    if (m % 64) {
        v.w[words - 1] &= (1ull << (m % 64)) - 1;
    }
    return v;
}

static std::vector<BitVec> random_rows(int n, int m, std::mt19937_64& rng) {
    std::vector<BitVec> rows;
    rows.reserve((size_t)n);
    for (int i = 0; i < n; ++i) {
        rows.push_back(random_bitvec(m, rng));
    }
    return rows;
}

static std::vector<BitVec> full_rank_rows(int m, std::mt19937_64& rng) {
    std::vector<BitVec> rows;
    rows.reserve((size_t)m);
    size_t words = (m + 63) / 64;
    for (int i = 0; i < m; ++i) {
        BitVec v = BitVec::make(m);
        for (size_t w = 0; w < words; ++w) v.w[w] = 0;
        int wi = i / 64;
        int bi = i & 63;
        v.w[wi] |= (1ull << bi);
        for (int j = i + 1; j < m; ++j) {
            if (rng() & 1ull) {
                int wj = j / 64;
                int bj = j & 63;
                v.w[wj] ^= (1ull << bj);
            }
        }
        rows.push_back(v);
    }
    return rows;
}

static std::vector<uint8_t> bitvec_to_bits(const BitVec& v, int m) {
    std::vector<uint8_t> out((size_t)m);
    for (int i = 0; i < m; ++i) {
        out[i] = (uint8_t)((v.w[i / 64] >> (i & 63)) & 1ull);
    }
    return out;
}

static bool gauss_mod2_solve(const std::vector<BitVec>& rows,
                             const std::vector<uint8_t>& rhs,
                             std::vector<uint8_t>& sol, int m) {
    int n = (int)rows.size();
    if (n == 0 || m <= 0) return false;

    size_t wlen = (m + 63) / 64;
    std::vector<std::vector<uint64_t>> M(n, std::vector<uint64_t>(wlen));
    std::vector<uint8_t> b = rhs;

    for (int i = 0; i < n; ++i) {
        for (size_t w = 0; w < wlen; ++w) {
            M[i][w] = rows[i].w[w];
        }
    }

    std::vector<int> pivot_row(m, -1);
    int row = 0;

    for (int col = 0; col < m && row < n; ++col) {
        size_t w = (size_t)col / 64;
        uint64_t mask = 1ull << (col & 63);

        int sel = -1;
        for (int i = row; i < n; ++i) {
            if (M[i][w] & mask) { sel = i; break; }
        }
        if (sel == -1) continue;

        if (sel != row) {
            std::swap(M[sel], M[row]);
            std::swap(b[sel], b[row]);
        }

        pivot_row[col] = row;

        for (int i = 0; i < n; ++i) {
            if (i == row) continue;
            if (M[i][w] & mask) {
                for (size_t j = 0; j < wlen; ++j) {
                    M[i][j] ^= M[row][j];
                }
                b[i] ^= b[row];
            }
        }
        row++;
    }

    for (int i = 0; i < n; ++i) {
        bool zero = true;
        for (size_t j = 0; j < wlen; ++j) {
            if (M[i][j] != 0) { zero = false; break; }
        }
        if (zero && b[i]) return false;
    }

    sol.assign((size_t)m, 0);
    for (int col = m - 1; col >= 0; --col) {
        int r = pivot_row[col];
        if (r < 0) { sol[col] = 0; continue; }

        uint8_t val = b[r];
        for (int c = col + 1; c < m; ++c) {
            if (M[r][(size_t)c / 64] & (1ull << (c & 63))) {
                val ^= sol[c];
            }
        }
        sol[col] = val;
    }
    return true;
}

int main() {
    std::cout << "- lpn test -\n";

    std::mt19937_64 rng(0x123456789abcdef0ull);

    {
        constexpr int m0 = 512;
        std::cout << "lp (no noise): m = " << m0 << "\n";

        BitVec secret0 = random_bitvec(m0, rng);
        auto rows0 = full_rank_rows(m0, rng);

        std::vector<uint8_t> y0(m0);
        for (int i = 0; i < m0; ++i) {
            y0[i] = bitvec_dot2(rows0[i], secret0);
        }

        std::vector<uint8_t> sol0;
        bool ok0 = gauss_mod2_solve(rows0, y0, sol0, m0);
        assert(ok0);

        auto s_bits = bitvec_to_bits(secret0, m0);
        for (int i = 0; i < m0; ++i) {
            assert(sol0[i] == s_bits[i]);
        }
        std::cout << "lp rec: ok\n";
    }


    // good params tbh
    constexpr int m = 4096;
    constexpr int n = 16384;
    constexpr double tau = 0.125;
    ///



    std::cout << "lpn: m = " << m << " n = " << n << " tau = " << tau << "\n";

    BitVec secret = random_bitvec(m, rng);
    auto rows = random_rows(n, m, rng);
    std::bernoulli_distribution noise(tau);

    std::vector<uint8_t> y(n);
    int wt_e = 0;
    int wt_y = 0;

    for (int i = 0; i < n; ++i) {
        uint8_t e = noise(rng) ? 1 : 0;
        uint8_t v = bitvec_dot2(rows[i], secret);
        uint8_t yi = v ^ e;
        y[i] = yi;
        wt_e += e;
        wt_y += yi;
    }

    double exp_e = tau * n;
    double sig_e = std::sqrt(n * tau * (1.0 - tau));
    double z_e = (wt_e - exp_e) / sig_e;
    assert(std::fabs(z_e) < 6.0);
    std::cout << "noise: wt = " << wt_e << " exp = " << exp_e << " z = " << z_e << "\n";

    double exp_y = 0.5 * n;
    double sig_y = std::sqrt(0.25 * n);
    double z_y = (wt_y - exp_y) / sig_y;
    assert(std::fabs(z_y) < 6.0);
    std::cout << "y: wt = " << wt_y << " exp = " << exp_y << " z = " << z_y << "\n";

    std::vector<uint8_t> sol;
    bool ok = gauss_mod2_solve(rows, y, sol, m);
    assert(!ok);

    std::cout << "noisy inc: ok\n";
    std::cout << "PASS\n";

    return 0;
}