# test_canon_tables.py -- record-level check against the canonical tables
# written by controlledburn-rs (tests/canon.rs): every run, edge, line and
# point record of both modes, fractions and lengths compared as f32 bits.
# The C++ core reads the same file (cpp/tests/test_canon_tables.cpp).

import pathlib
import struct

import pytest

from controlledburn import _core

CANON = pathlib.Path(__file__).resolve().parents[2] / "fixtures" / "canon-tables.txt"


def _f64(h):
    return struct.unpack(">d", bytes.fromhex(h))[0]


def _f32_hex(v):
    return struct.pack(">f", float(v)).hex()


def load_cases():
    cases, cur, mode = [], None, None
    for line in CANON.read_text().splitlines():
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "case":
            cur = {"name": t[1], "wkb": [], "tables": {}}
        elif t[0] == "grid":
            cur["grid"] = tuple(_f64(h) for h in t[1:5]) + (int(t[5]), int(t[6]))
        elif t[0] == "wkb":
            cur["wkb"].append(bytes.fromhex(t[1]))
        elif t[0] == "mode":
            mode = t[1]
            cur["tables"][mode] = []
        elif t[0] == "end":
            cases.append(cur)
        else:
            cur["tables"][mode].append(line.strip())
    return cases


def records(raw):
    out = []
    r = raw["runs"]
    out += [f"run {a} {b} {c} {d}" for a, b, c, d in
            zip(r["row"], r["col_start"], r["col_end"], r["id"])]
    e = raw["edges"]
    out += [f"edge {a} {b} {_f32_hex(c)} {d}" for a, b, c, d in
            zip(e["row"], e["col"], e["fraction"], e["id"])]
    li = raw["lines"]
    out += [f"line {a} {b} {_f32_hex(c)} {d}" for a, b, c, d in
            zip(li["row"], li["col"], li["length"], li["id"])]
    p = raw["points"]
    out += [f"point {a} {b} {c}" for a, b, c in zip(p["row"], p["col"], p["id"])]
    return out


if not CANON.exists():
    pytest.skip("fixtures/ directory not found", allow_module_level=True)

CASES = load_cases()


@pytest.mark.parametrize("mode", ["coverage", "approx"])
def test_canonical_tables(mode):
    bad = []
    for c in CASES:
        # the canon is parity mode: no run coalescing, no part merging
        got = records(_core.burn_wkb(c["wkb"], *c["grid"], mode,
                                     coalesce=False, merge_parts=False))
        if got != c["tables"][mode]:
            bad.append(c["name"])
    assert len(CASES) == 240
    assert not bad, f"{len(bad)} cases differ from canon-tables.txt: {bad[:5]}"
