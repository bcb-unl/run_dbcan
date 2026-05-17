import numpy as np

from dbcan.plot.expression_plots import build_sem_trap_trace, gene_bp_to_pixel


def test_gene_bp_to_pixel_spans():
    starts = [100, 300]
    ends = [200, 400]
    x0, x1, shift, pix = gene_bp_to_pixel(starts, ends, width=1000)
    assert shift == 100
    assert x0[0] == 0
    assert abs(x1[0] - 1000 * (200 - 100) / (400 - 100)) < 1e-6
    assert x1[1] == 1000


def test_step_trace_intergenic_baseline():
    xs, ys = build_sem_trap_trace([0, 200], [100, 300], [5.0, 10.0], baseline=0.0)
    assert xs[0] == 0 and xs[1] == 100
    assert ys[0] == 5.0 and ys[1] == 5.0
    assert ys[2] == 0.0 and ys[3] == 0.0
    assert xs[4] == 200
    assert ys[-1] == 10.0


def test_step_trace_single_gene():
    xs, ys = build_sem_trap_trace([10], [50], [3.0])
    assert len(xs) == 2
    assert list(ys) == [3.0, 3.0]
