"""
Check the outputs of test_gaussfilt.f90 (gaussfilt_out/): gauss_filt with the default method 'rsf' on short
arrays and for small to huge sigma, in 1D, 2D and 3D, in single and double precision.

Along every axis, the reference is the direct convolution of the input, extended by repeated mirror
reflection (period 2(n - 1), the edge sample not repeated), with the kernel that gauss_filt uses:
  sigma >= 1       the two-term recursive spectral kernel (Sugimoto and Kamata, 2013),
                   h(u) = [1 + s1 cos(phi u) + s2 cos(2 phi u)]/(2r + 1) for |u| <= r, with phi = 2 pi/(2r + 1),
                   s_k = 2 exp(-sigma^2 phi^2 k^2/2), r = nint(3 sigma + 0.5) (nint(3 sigma + 0.3) for sigma < 4)
  0 < sigma < 1    the normalized sampled Gaussian for |u| <= ceil(6 sigma)
The kernel is folded onto one period of the extension, which is exact for any radius. Two exact conditions
are checked as well: the output equals the input for sigma < 1/15 (to 1e-40 of the largest input value: the
neighbor weight of the sampled Gaussian, at most 2e-49, is representable in double precision), and a
constant input stays constant.

Run: python check_gaussfilt.py [directory, default gaussfilt_out]
Prints the largest errors (relative to the largest input value) by sigma range and writes gaussfilt.png
(1D) and gaussfilt_2d3d.png (2D and 3D examples).
"""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# On short arrays the recursion runs over about 6 sigma samples, and its rounding grows with that length
# (about 2e-11 in double precision at sigma = 1e6)
TOL = {'float': 1.0e-5, 'double': 1.0e-10}
DTYPE = {'float': np.float32, 'double': np.float64}


def fnint(x):
    """Fortran nint: round half away from zero."""
    return int(np.sign(x) * np.floor(abs(x) + 0.5))


def kernel(sigma, prec):
    """Offsets and weights of the kernel of gauss_filt; the radius is computed in the precision of the call."""
    t = DTYPE[prec]
    s = t(sigma)
    if s >= 1.0:
        r = fnint(t(3.0) * s - t(0.2) + t(0.5)) if s < 4.0 else fnint(t(3.0) * s + t(0.5))
        phi = 2.0 * np.pi / (2 * r + 1)
        s = float(s)
        u = np.arange(-r, r + 1)
        return u, (1.0 + 2.0 * np.exp(-0.5 * s ** 2 * phi ** 2) * np.cos(phi * u)
                   + 2.0 * np.exp(-2.0 * s ** 2 * phi ** 2) * np.cos(2.0 * phi * u)) / (2 * r + 1)
    r = int(np.ceil(t(6.0) * s))
    u = np.arange(-r, r + 1)
    with np.errstate(under='ignore'):
        g = np.exp(-0.5 * u ** 2 / float(s) ** 2)
    return u, g / g.sum()


def filt_axis(w, sigma, prec, axis):
    n = w.shape[axis]
    if sigma <= 0 or n == 1:
        return w
    u, h = kernel(sigma, prec)
    p = 2 * (n - 1)
    folded = np.bincount(np.mod(u, p), weights=h, minlength=p)
    v = np.moveaxis(w, axis, 0)
    m = np.arange(p)
    ext = v[(n - 1) - np.abs((n - 1) - m)]
    idx = (np.arange(n)[:, None] + m[None, :]) % p
    y = np.tensordot(ext[idx], folded, axes=([1], [0]))
    return np.moveaxis(y, 0, axis)


def reference(w, sigmas, prec):
    for ax, s in enumerate(sigmas):
        w = filt_axis(w, s, prec, ax)
    return w


def read(path, shape, dt):
    return np.fromfile(path, dtype=dt).astype(np.float64).reshape(shape[::-1]).T


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else 'gaussfilt_out'
    rows = []
    for line in open(d + '/cases.txt'):
        f = line.split()
        k, nd = int(f[0]), int(f[1])
        shape = tuple(int(v) for v in f[2:2 + nd])
        sig = [float(np.float32(v)) for v in f[5:5 + nd]]
        w = read('%s/case_%d_in.bin' % (d, k), shape, np.float32)
        scale = max(np.abs(w).max(), 1.0e-30)
        q = dict(k=k, nd=nd, shape=shape, sig=sig, signal=f[8], w=w)
        for prec, suf in (('float', 'out'), ('double', 'outd')):
            y = read('%s/case_%d_%s.bin' % (d, k, suf), shape, DTYPE[prec])
            q[prec] = y
            q['err_' + prec] = np.abs(y - reference(w, sig, prec)).max() / scale
            q['finite_' + prec] = bool(np.all(np.isfinite(y)))
        active = [s for s, n in zip(sig, shape) if s > 0 and n > 1]
        q['identity'] = all(s < 1.0 / 15.0 for s in active)
        q['range'] = ('none' if not active else
                      'sigma < 1/15' if max(active) < 1.0 / 15.0 else
                      '1/15 <= sigma < 1' if max(active) < 1.0 else
                      '1 <= sigma <= 100' if max(active) <= 100.0 else 'sigma >= 1000')
        rows.append(q)

    print('%d cases in %s; errors are max |gauss_filt - reference| / max |input|' % (len(rows), d))
    print('%-22s %-6s %6s   %-14s %-14s' % ('largest sigma', 'dim', 'cases', 'float', 'double'))
    failed = []
    for rng in ('none', 'sigma < 1/15', '1/15 <= sigma < 1', '1 <= sigma <= 100', 'sigma >= 1000'):
        for nd in (1, 2, 3):
            c = [q for q in rows if q['range'] == rng and q['nd'] == nd]
            if not c:
                continue
            ef = max(q['err_float'] for q in c)
            ed = max(q['err_double'] for q in c)
            print('%-22s %-6s %6d   %-14s %-14s' % (rng, '%dD' % nd, len(c),
                                                     '%.1e %s' % (ef, 'pass' if ef <= TOL['float'] else 'FAIL'),
                                                     '%.1e %s' % (ed, 'pass' if ed <= TOL['double'] else 'FAIL')))
            failed += [q for q in c if q['err_float'] > TOL['float'] or q['err_double'] > TOL['double']]

    ident = [q for q in rows if q['identity']]
    exact = [q for q in ident if all(np.abs(q[p] - q['w']).max() <= 1.0e-40 * max(np.abs(q['w']).max(), 1.0e-30)
                                     for p in ('float', 'double'))]
    bitwise = [q for q in ident if np.array_equal(q['float'], q['w'])]
    print('identity (all filtered axes with sigma < 1/15): output == input to 1e-40 in %d of %d cases '
          '(bit-identical in single precision in %d)' % (len(exact), len(ident), len(bitwise)))
    const = [q for q in rows if q['signal'] == 'const']
    print('constant input: max |output - 1| = %.1e (float), %.1e (double), %d cases' % (
        max(np.abs(q['float'] - 1.0).max() for q in const), max(np.abs(q['double'] - 1.0).max() for q in const), len(const)))
    finite = all(q['finite_float'] and q['finite_double'] for q in rows)
    print('all outputs finite: %s' % finite)
    if failed:
        print('FAILED cases (tolerance float %.0e, double %.0e):' % (TOL['float'], TOL['double']))
        for q in failed[:20]:
            print('  case %d: %dD %s sigma %s %s: float %.1e, double %.1e' % (q['k'], q['nd'], q['shape'], q['sig'],
                                                                              q['signal'], q['err_float'], q['err_double']))
    ok = not failed and len(exact) == len(ident) and finite
    print('RESULT: %s' % ('PASS' if ok else 'FAIL'))

    plot(rows, 'gaussfilt.png')
    plot_2d3d(rows, 'gaussfilt_2d3d.png')
    return 0 if ok else 1


def plot(rows, path):
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 14})
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.2))
    one = [q for q in rows if q['nd'] == 1 and q['shape'][0] > 1]
    sig = sorted(set(q['sig'][0] for q in one))
    for prec, col in (('float', 'red'), ('double', 'blue')):
        e = [max(max(q['err_' + prec] for q in one if q['sig'][0] == s), 1.0e-17) for s in sig]
        axs[0].loglog(sig, e, 'o-', color=col, lw=2, ms=6, label=prec.capitalize())
        axs[0].axhline(TOL[prec], color=col, lw=1, ls='--')
    for x in (1.0 / 15.0, 1.0):
        axs[0].axvline(x, color='black', lw=1, ls=':')
    axs[0].set_xlabel('$\\sigma$')
    axs[0].set_ylabel('Max error / max |input|')
    axs[0].set_title('(a) 1D, n = 2 to 16, all inputs', loc='left')
    axs[0].legend(frameon=False)
    q10 = [q for q in one if q['shape'][0] == 10 and q['signal'] == 'random']
    x = np.arange(1, 11)
    axs[1].plot(x, q10[0]['w'], 'o', color='gray', ms=8, label='Input')
    for s, col in zip((0.25, 1.0, 10.0, 1.0e3, 1.0e6), ('black', 'blue', 'red', 'green', 'magenta')):
        q = [r for r in q10 if np.isclose(r['sig'][0], s)][0]
        axs[1].plot(x, q['float'], '-', color=col, lw=2, label='$\\sigma$ = %g' % s)
        axs[1].plot(x, reference(q['w'], q['sig'], 'float'), 's', color=col, ms=4, mfc='none')
    axs[1].set_xlabel('Sample')
    axs[1].set_ylabel('Value')
    axs[1].set_title('(b) n = 10: gauss_filt (lines), reference (squares)', loc='left')
    axs[1].legend(frameon=False, fontsize=12, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    print('wrote %s' % path)


def plot_2d3d(rows, path):
    """Input, gauss_filt (single precision), reference and difference for a few small 2D and 3D arrays."""
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 13})
    sel = [(2, (10, 10), (10.0, 10.0), 'spike1'), (2, (7, 12), (10.0, 3.0), 'random'),
           (2, (3, 20), (0.25, 3.0), 'random'), (3, (10, 10, 10), (10.0, 10.0, 10.0), 'spike1'),
           (3, (4, 5, 6), (0.05, 0.25, 0.5), 'random'), (3, (6, 8, 10), (1.0e5, 1.0e5, 1.0e5), 'random')]
    fig, axs = plt.subplots(len(sel), 4, figsize=(15, 3.1 * len(sel)))
    for i, (nd, shape, sig, kind) in enumerate(sel):
        q = [r for r in rows if r['nd'] == nd and r['shape'] == shape and r['signal'] == kind
             and np.allclose(r['sig'], sig, rtol=1.0e-6)][0]
        cut = (lambda z: z) if nd == 2 else (lambda z: z[:, :, 0])
        w, v, a = cut(q['w']), cut(q['float']), cut(reference(q['w'], q['sig'], 'float'))
        lim = (min(v.min(), a.min()), max(v.max(), a.max()))
        if lim[1] - lim[0] < 1.0e-12:
            lim = (lim[0] - 1.0e-6, lim[1] + 1.0e-6)
        for j, (z, t) in enumerate([(w, 'Input'), (v, 'gauss_filt'), (a, 'Reference'), (v - a, 'gauss_filt - reference')]):
            ax = axs[i, j]
            if j == 3:
                c = max(np.abs(z).max(), 1.0e-12)
                im = ax.imshow(z, cmap='bwr', vmin=-c, vmax=c, aspect='auto', interpolation='none')
                fmt = '%.1e'
            else:
                vmin, vmax = (w.min(), w.max()) if j == 0 else lim
                im = ax.imshow(z, cmap='jet', vmin=vmin, vmax=vmax, aspect='auto', interpolation='none')
                fmt = '%.2f' if j == 0 else None
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, format=fmt)
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_title(t)
        axs[i, 0].set_ylabel('%s %s\n$\\sigma$ = %s%s' % ('x'.join(str(n) for n in shape), kind,
                             ', '.join('%g' % s for s in sig), '' if nd == 2 else '\n(slice 1 of axis 3)'), fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    print('wrote %s' % path)


if __name__ == '__main__':
    sys.exit(main())
