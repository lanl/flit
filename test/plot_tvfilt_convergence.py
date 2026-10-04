#!/usr/bin/env python3
"""
Plot the convergence curves written by test_tvfilt.f90.

Usage:
    python plot_tvfilt_convergence.py [--input test_tvfilt_convergence.txt]
                                      [--out test_tvfilt_convergence.png]
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

SURFACE, INK, INK2, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#e1e0d9', '#c3c2b7'
COLORS = {('tv_filt', 0.0): '#2a78d6', ('tgpv_filt', 1.0): '#eb6834', ('tgpv_filt', 0.5): '#1baf7a'}
LABELS = {('tv_filt', 0.0): 'TV', ('tgpv_filt', 1.0): 'TGpV, p = 1', ('tgpv_filt', 0.5): 'TGpV, p = 0.5'}
TOL = {'constant': 1e-12, 'noisy': 0.6}
FLOOR = 1e-17


def read_curves(path):
    iterations = None
    curves = {}
    with open(path) as f:
        for line in f:
            if line.startswith('# iterations'):
                iterations = np.array([int(v) for v in line.split()[2:]])
            elif line.strip() and not line.startswith('#'):
                v = line.split()
                curves[(v[0], int(v[1]), v[2], float(v[3]))] = np.array([float(e) for e in v[4:]])
    return iterations, curves


def main():
    parser = argparse.ArgumentParser(description='Plot the test_tvfilt convergence curves.')
    parser.add_argument('--input', default='test_tvfilt_convergence.txt')
    parser.add_argument('--out', default='test_tvfilt_convergence.png')
    args = parser.parse_args()

    if any(f.name == 'Arial' for f in font_manager.fontManager.ttflist):
        plt.rcParams['font.family'] = 'Arial'
    plt.rcParams.update({
        'font.size': 11, 'axes.titlesize': 11, 'legend.fontsize': 10,
        'xtick.labelsize': 10, 'ytick.labelsize': 10,
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
        'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': AXIS,
        'xtick.color': INK2, 'ytick.color': INK2,
    })

    iterations, curves = read_curves(args.input)
    fig, axs = plt.subplots(2, 3, figsize=(7.5, 5.0), sharex=True, layout='constrained')
    for row, image in enumerate(['constant', 'noisy']):
        for col, nd in enumerate([1, 2, 3]):
            ax = axs[row, col]
            ax.grid(True, color=GRID, linewidth=0.8)
            ax.set_axisbelow(True)
            for side in ['top', 'right']:
                ax.spines[side].set_visible(False)
            ax.axhline(TOL[image], color=AXIS, linewidth=1.0)
            for key, color in COLORS.items():
                e = curves[(image, nd, key[0], key[1])]
                if image == 'constant':
                    e = np.maximum(e, FLOOR)
                ax.plot(iterations, e, color=color, lw=2, solid_capstyle='round', label=LABELS[key])
                ax.plot(iterations, e, 'o', ms=5, color=color, mec=SURFACE, mew=1.2)
            ax.set_xscale('log')
            if image == 'constant':
                ax.set_yscale('log')
                ax.set_ylim(FLOOR/3, 1e-1)
                ax.set_title(f'({"abc"[col]})  Constant image, {nd}D', loc='left')
            else:
                ax.set_ylim(0, 1)
                ax.set_title(f'({"def"[col]})  Noisy image, {nd}D', loc='left')
                ax.set_xlabel('Iterations')
            if col == 0:
                ax.set_ylabel('max|u − c| / c' if image == 'constant' else 'rms(u − truth) / rms(g − truth)')
    axs[0, 0].legend(loc='center right', frameon=False)
    fig.savefig(args.out, dpi=200)
    print(f'Wrote {args.out}')


if __name__ == '__main__':
    main()
