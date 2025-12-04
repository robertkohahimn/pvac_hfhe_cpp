#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cstring>
#include <pvac/pvac.hpp>

using namespace pvac;

static bool test_sha256_abc() {
    const char *msg = "abc";
    uint8_t out[32];

    sha256_bytes(msg, std::strlen(msg), out);

    const uint8_t ref[32] = {
        0xba,0x78,0x16,0xbf,0x8f,0x01,0xcf,0xea,
        0x41,0x41,0x40,0xde,0x5d,0xae,0x22,0x23,
        0xb0,0x03,0x61,0xa3,0x96,0x17,0x7a,0x9c,
        0xb4,0x10,0xff,0x61,0xf2,0x00,0x15,0xad
    };

    return std::memcmp(out, ref, 32) == 0;
}

static bool test_xof_basic() {
    std::vector<uint64_t> seed = {1,2,3,4};
    XofShake x1, x2;

    x1.init(std::string("test"), seed);
    x2.init(std::string("test"), seed);

    uint8_t a[64], b[64];

    for (int i=0;i<4;i++) {
        uint64_t w1 = x1.take_u64();
        uint64_t w2 = x2.take_u64();

        if (w1 != w2) return false;

        std::memcpy(a + 8*i, &w1, 8);
        std::memcpy(b + 8*i, &w2, 8);
    }
    XofShake x3;
    x3.init(std::string("test2"), seed);




    uint64_t w3 = x3.take_u64();
    return (w3 != *reinterpret_cast<uint64_t*>(a));
}

static int hamming_64(uint64_t x) {
    int c = 0;
    while (x) { x &= x-1; c++; }
    return c;
}

static bool test_prf_R_domains() {
    Params prm;
    PubKey pk;
    SecKey sk;
    keygen(prm, pk, sk);
    RSeed seed;
    seed.ztag = csprng_u64();
    seed.nonce = make_nonce128();



    Fp r1 = prf_R_core(pk, sk, seed, Dom::PRF_R1);
    Fp r2 = prf_R_core(pk, sk, seed, Dom::PRF_R2);
    if (r1.lo == r2.lo && r1.hi == r2.hi) return false;
    int hw_lo = hamming_64(r1.lo ^ r2.lo);
    int hw_hi = hamming_64(r1.hi ^ r2.hi);
    int hw = hw_lo + hw_hi;
    return hw > 40 && hw < 88;
}

int main() {
    bool ok1 = test_sha256_abc();
    bool ok2 = test_xof_basic();
    bool ok3 = test_prf_R_domains();
    std::cout << "- prf/hash tests -\n";
    std::cout << "sha256(abc): " << (ok1 ? "ok" : "FAIL") << "\n";
    std::cout << "xof: " << (ok2 ? "ok" : "FAIL") << "\n";
    std::cout << "prf_R domains: " << (ok3 ? "ok" : "FAIL") << "\n";
    bool all = ok1 && ok2 && ok3;
    std::cout << "\nresult: " << (all ? "PASS" : "FAIL") << "\n";
    return all ? 0 : 1;
}