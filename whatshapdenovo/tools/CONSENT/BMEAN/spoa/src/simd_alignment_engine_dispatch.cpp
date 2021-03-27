/*!
 * @file simd_alignment_engine_dispatch.cpp
 *
 * @brief Instantiation of different SIMD engines
 */

 #include "simd_alignment_engine_impl.hpp"

 #if defined(__AVX2__)
 #define ARCH Arch::avx2
 #elif defined (__SSE4_1__)
 #define ARCH Arch::sse4_1
 #else
 #define ARCH Arch::sse2
 #endif


namespace spoa{

template class SimdAlignmentEngine<ARCH>;

template
std::unique_ptr<AlignmentEngine> createSimdAlignmentEngine<ARCH>(AlignmentType type,
    AlignmentSubtype subtype, std::int8_t m, std::int8_t n, std::int8_t g,
    std::int8_t e, std::int8_t q, std::int8_t c);

}

