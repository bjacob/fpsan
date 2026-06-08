# fpsan

A personal repository for experiments, study, and documentation around
**FPSan**, the floating-point sanitizer invented in the
[Triton project](https://github.com/triton-lang/triton). Nothing here is
original: FPSan and the ideas explored are Triton's. See the authors' blog post,
["Schanuel's conjecture and the semantics of FPSan"](https://cp4space.hatsya.com/2026/05/03/schanuels-conjecture-and-the-semantics-of-fpsan/),
and the [Triton FPSan docs](https://github.com/triton-lang/triton/blob/main/docs/programming-guide/chapter-3/fpsan.rst).

The real, complete implementation lives at
[ROCm/hip-fpsan](https://github.com/ROCm/hip-fpsan).

- [understanding-fpsan.md](understanding-fpsan.md) — the construction derived from first principles.
- [narrow-floats.md](narrow-floats.md) — the same recipe worked by hand on 4- and 8-bit floats.
- `fpsan_f32.hpp`, `fpsan_generic.hpp`, `demo.cpp`, `narrow_demo.cpp` — the standalone toy code.
