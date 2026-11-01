# _core.py -- the binding surface used by the package API.
#
# Forwards to the native module built from the controlledburn Rust crate
# (rust/src/lib.rs), coercing array arguments to the contiguous dtypes the
# native functions take. Same names and signatures as the former pybind11
# module over the C++ core, so the package API above it is unchanged.
#
# Copyright (c) 2025 Michael Sumner
# Licensed under Apache License 2.0

import numpy as np

from . import _native

__doc__ = _native.__doc__


def _i64(a):
    return np.ascontiguousarray(a, dtype=np.int64)


def _mask(valid):
    return None if valid is None else np.ascontiguousarray(valid, dtype=np.uint8)


def burn_wkb(wkb_list, xmin, ymin, xmax, ymax, ncol, nrow, mode="coverage"):
    return _native.burn_wkb(wkb_list, float(xmin), float(ymin), float(xmax),
                            float(ymax), int(ncol), int(nrow), mode)


def burn_wkb_arrow(values, offsets, valid, xmin, ymin, xmax, ymax, ncol, nrow,
                   mode="coverage"):
    return _native.burn_wkb_arrow(values, _i64(offsets), _mask(valid),
                                  float(xmin), float(ymin), float(xmax),
                                  float(ymax), int(ncol), int(nrow), mode)


def bbox_wkb(wkb_list):
    return _native.bbox_wkb(wkb_list)


def bbox_wkb_arrow(values, offsets, valid):
    return _native.bbox_wkb_arrow(values, _i64(offsets), _mask(valid))


def materialize(runs_row, runs_col_start, runs_col_end, runs_id,
                edges_row, edges_col, edges_fraction, edges_id,
                ncol, nrow, values, fn, edge_policy, threshold, background):
    i4 = lambda a: np.ascontiguousarray(a, dtype=np.int32)  # noqa: E731
    return _native.materialize(
        i4(runs_row), i4(runs_col_start), i4(runs_col_end), i4(runs_id),
        i4(edges_row), i4(edges_col),
        np.ascontiguousarray(edges_fraction, dtype=np.float32), i4(edges_id),
        int(ncol), int(nrow), [float(v) for v in values], fn, edge_policy,
        float(threshold), float(background))
