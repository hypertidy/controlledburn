// no_fma.h -- turn off floating-point contraction (fused multiply-add)
// for the translation unit that includes it. Include it before any other
// header in every core source file.
//
// The Rust crate controlledburn-rs, the canonical implementation, never
// fuses a multiply and an add. clang (arm64 by default) and GCC (any
// target with FMA, e.g. aarch64 or -march=x86-64-v3) otherwise do, which
// moves run boundaries and line records on knife-edge inputs, so results
// would depend on the platform. cpp/tests/test_canon_tables.cpp checks
// this against fixtures/canon-tables.txt. Every core translation unit
// must agree: one contracted file feeding an uncontracted one is worse
// than either alone.
//
// This is a pragma rather than -ffp-contract=off so that every build
// route (CMake, R's Makevars, Python) gets it without compiler flags.

#ifndef CONTROLLEDBURN_NO_FMA_H
#define CONTROLLEDBURN_NO_FMA_H

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#elif defined(_MSC_VER)
#pragma fp_contract(off)
#endif

#endif
