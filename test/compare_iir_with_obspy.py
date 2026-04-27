#!/usr/bin/env python3
"""
Compare FLIT IIR against an ObsPy IIR for five filter modes: 
lowpass, highpass, bandpass, bandstop, notch.

Assumed executable interface:

    ./exec type fs flow fhigh order zerophase

Examples:
    ./exec bandpass 200.0 10.0 20.0 4 1
    ./exec lowpass  200.0 0.0  30.0 4 1
    ./exec highpass 200.0 5.0  0.0  4 1
    ./exec bandstop 200.0 45.0 55.0 4 1
    ./exec notch    200.0 60.0 2.0  4 1

Assumptions:
    1. This script writes the input signal to x.txt by default.
    2. Your Fortran executable reads x.txt by default.
    3. Your Fortran executable writes y.txt by default.
    4. Text vector format is plain one-column data with no n_samples header:
           x1
           x2
           x3
           ...

Use --input-file and --output-file if your executable uses different names.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import iirfilter, sosfilt, chirp


def write_vec(path: str | Path, x: np.ndarray) -> None:
    """
    Write plain one-column signal file:

        x1
        x2
        x3
        ...

    No n_samples header.
    """
    path = Path(path)
    np.savetxt(path, np.asarray(x, dtype=np.float64), fmt="%.17e")


def read_vec(path: str | Path) -> np.ndarray:
    """
    Read plain one-column signal file:

        y1
        y2
        y3
        ...

    No n_samples header.
    """
    path = Path(path)
    y = np.loadtxt(path, dtype=np.float64)

    y = np.asarray(y, dtype=np.float64)
    if y.ndim == 0:
        y = y.reshape(1)

    return y


def colored_noise(n: int, alpha: float, rng: np.random.Generator) -> np.ndarray:
    x = rng.standard_normal(n)
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n)

    scale = np.ones_like(freqs)
    scale[0] = 0.0
    scale[1:] = 1.0 / np.maximum(freqs[1:], 1.0 / n) ** (0.5 * alpha)

    y = np.fft.irfft(X * scale, n=n)
    y -= np.mean(y)
    y /= max(np.std(y), 1.0e-30)
    return y


def make_random_signal(fs: float, duration: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)

    n = int(round(duration * fs))
    t = np.arange(n, dtype=np.float64) / fs

    x = np.zeros(n, dtype=np.float64)

    x += 0.50 * colored_noise(n, alpha=1.0, rng=rng)
    x += 0.12 * rng.standard_normal(n)
    x += 0.25 * colored_noise(n, alpha=2.0, rng=rng)

    for _ in range(30):
        center = rng.uniform(t[0], t[-1])
        width = rng.uniform(0.01, 0.12)
        amp = rng.normal(0.0, 0.8)
        x += amp * np.exp(-0.5 * ((t - center) / width) ** 2)

    for _ in range(5):
        center = rng.uniform(0.15 * duration, 0.85 * duration)
        length = rng.uniform(0.7, 2.3)
        f0 = rng.uniform(2.0, 12.0)
        f1 = min(rng.uniform(25.0, 80.0), 0.85 * 0.5 * fs)
        amp = rng.uniform(0.15, 0.45)

        tau = t - (center - 0.5 * length)
        local_t = np.clip(tau, 0.0, length)
        c = chirp(local_t, f0=f0, f1=f1, t1=length, method="linear")
        window = np.exp(-0.5 * ((t - center) / (0.22 * length)) ** 2)
        x += amp * c * window

    # Coherent components chosen to exercise the default filters for fs=200.
    x += 0.60 * np.sin(2.0 * np.pi * 6.0 * t + 0.3 * np.sin(2.0 * np.pi * 0.15 * t))
    x += 0.45 * np.sin(2.0 * np.pi * 18.0 * t + rng.uniform(0.0, 2.0 * np.pi))
    x += 0.30 * np.sin(2.0 * np.pi * min(70.0, 0.7 * 0.5 * fs) * t + rng.uniform(0.0, 2.0 * np.pi))

    # Strong narrow spike/tone for notch removal.
    notch_freq = min(60.0, 0.6 * 0.5 * fs)
    x += 1.10 * np.sin(2.0 * np.pi * notch_freq * t + rng.uniform(0.0, 2.0 * np.pi))

    x -= np.mean(x)
    x /= max(np.std(x), 1.0e-30)

    return t, x


def obspy_style_filter(x: np.ndarray, sos: np.ndarray, zerophase: bool) -> np.ndarray:
    y = sosfilt(sos, x)
    if zerophase:
        y = sosfilt(sos, y[::-1])[::-1]
    return y


def scipy_reference_sos(case: dict, fs: float, order: int, ftype: str = "butter") -> np.ndarray:
    nyq = 0.5 * fs
    kind = case["type"]
    flow = case["flow"]
    fhigh = case["fhigh"]

    if kind == "lowpass":
        Wn = fhigh / nyq
        btype = "lowpass"
    elif kind == "highpass":
        Wn = flow / nyq
        btype = "highpass"
    elif kind == "bandpass":
        Wn = [flow / nyq, fhigh / nyq]
        btype = "bandpass"
    elif kind == "bandstop":
        Wn = [flow / nyq, fhigh / nyq]
        btype = "bandstop"
    elif kind == "notch":
        center = flow
        width = fhigh
        f1 = center - 0.5 * width
        f2 = center + 0.5 * width
        Wn = [f1 / nyq, f2 / nyq]
        btype = "bandstop"
    else:
        raise ValueError(f"Unknown filter type: {kind}")

    return iirfilter(N=order, Wn=Wn, btype=btype, ftype=ftype, output="sos")


def run_fortran_exec(
    exe: str | Path,
    case: dict,
    fs: float,
    order: int,
    zerophase: bool,
    output_file: str | Path,
) -> np.ndarray:
    exe = Path(exe).expanduser().resolve()
    if not exe.exists():
        raise FileNotFoundError(f"Cannot find executable: {exe}")

    zflag = "1" if zerophase else "0"

    cmd = [
        str(exe),
        str(case["type"]),
        f"{fs:.17g}",
        f"{case['flow']:.17g}",
        f"{case['fhigh']:.17g}",
        str(order),
        zflag,
    ]

    print("Running:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)

    return read_vec(output_file)


def amplitude_spectrum(x: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    n = x.size
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    a = np.abs(np.fft.rfft(x)) / n
    return f, a


def plot_all(
    t: np.ndarray,
    x: np.ndarray,
    results: list[dict],
    fs: float,
    out: str | Path,
    show: bool,
) -> None:
    out = Path(out)
    ncase = len(results)

    fig, axes = plt.subplots(
        nrows=1 + 2 * ncase,
        ncols=1,
        figsize=(14, 3.0 + 2.2 * ncase),
        sharex=True,
    )

    axes[0].plot(t, x, linewidth=0.8)
    axes[0].set_title("Input random/nonstationary signal")
    axes[0].set_ylabel("amp.")
    axes[0].grid(True, alpha=0.3)

    for i, r in enumerate(results):
        ax_y = axes[1 + 2 * i]
        ax_e = axes[1 + 2 * i + 1]

        ax_y.plot(t, r["y_ref"], linewidth=1.2, label="ObsPy IIR")
        ax_y.plot(t, r["y_fortran"], "--", linewidth=0.9, label="FLIT IIR")
        ax_y.set_title(
            f"{r['label']} | max abs={r['max_abs']:.3e}, rel L2={r['rel_l2']:.3e}",
            fontsize=10,
        )
        ax_y.set_ylabel("amp.")
        ax_y.grid(True, alpha=0.3)

        if i == 0:
            ax_y.legend(loc="upper right")

        ax_e.plot(t, r["err"], linewidth=0.8)
        ax_e.set_ylabel("error")
        ax_e.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time [s]")
    fig.suptitle("FLIT vs. ObsPy IIR filtering", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(out, dpi=180)
    print(f"Saved time-domain comparison: {out}")

    out_spectrum = out.with_name(out.stem + "_spectra.png")
    fig2, axes2 = plt.subplots(
        nrows=1 + ncase,
        ncols=1,
        figsize=(14, 3.0 + 2.0 * ncase),
        sharex=True,
    )

    f, a = amplitude_spectrum(x, fs)
    axes2[0].semilogy(f, a + 1.0e-16, linewidth=1.0)
    axes2[0].set_title("Input amplitude spectrum")
    axes2[0].set_ylabel("amp.")
    axes2[0].grid(True, alpha=0.3)

    for i, r in enumerate(results):
        ax = axes2[i + 1]
        f, a_ref = amplitude_spectrum(r["y_ref"], fs)
        _, a_fortran = amplitude_spectrum(r["y_fortran"], fs)

        ax.semilogy(f, a_ref + 1.0e-16, linewidth=1.2, label="ObsPy IIR")
        ax.semilogy(f, a_fortran + 1.0e-16, "--", linewidth=0.9, label="FLIT IIR")
        ax.set_title(r["label"], fontsize=10)
        ax.set_ylabel("amp.")
        ax.grid(True, alpha=0.3)

        if i == 0:
            ax.legend(loc="upper right")

    axes2[-1].set_xlabel("Frequency [Hz]")
    fig2.suptitle("Amplitude spectra comparison", fontsize=14)
    fig2.tight_layout(rect=[0, 0, 1, 0.975])
    fig2.savefig(out_spectrum, dpi=180)
    print(f"Saved spectrum comparison: {out_spectrum}")

    if show:
        plt.show()


def default_cases(fs: float) -> list[dict]:
    nyq = 0.5 * fs

    if nyq <= 80.0:
        return [
            {"type": "lowpass",  "flow": 0.0,       "fhigh": 0.30 * nyq, "label": f"lowpass < {0.30*nyq:.3g} Hz"},
            {"type": "highpass", "flow": 0.08*nyq,  "fhigh": 0.0,        "label": f"highpass > {0.08*nyq:.3g} Hz"},
            {"type": "bandpass", "flow": 0.10*nyq,  "fhigh": 0.35*nyq,   "label": f"bandpass {0.10*nyq:.3g}-{0.35*nyq:.3g} Hz"},
            {"type": "bandstop", "flow": 0.35*nyq,  "fhigh": 0.45*nyq,   "label": f"bandstop {0.35*nyq:.3g}-{0.45*nyq:.3g} Hz"},
            {"type": "notch",    "flow": 0.60*nyq,  "fhigh": 0.03*nyq,   "label": f"notch center {0.60*nyq:.3g} Hz, width {0.03*nyq:.3g} Hz"},
        ]

    return [
        {"type": "lowpass",  "flow": 0.0,  "fhigh": 30.0, "label": "lowpass < 30 Hz"},
        {"type": "highpass", "flow": 5.0,  "fhigh": 0.0,  "label": "highpass > 5 Hz"},
        {"type": "bandpass", "flow": 10.0, "fhigh": 20.0, "label": "bandpass 10-20 Hz"},
        {"type": "bandstop", "flow": 45.0, "fhigh": 55.0, "label": "bandstop 45-55 Hz"},
        {"type": "notch",    "flow": 60.0, "fhigh": 2.0,  "label": "notch center 60 Hz, width 2 Hz"},
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Automatically compare FLIT IIR with ObsPy IIR. "
    )

    parser.add_argument("--exe", default="./exec", help="Path to your FLIT executable.")
    parser.add_argument("--fs", type=float, default=200.0, help="Sampling frequency in Hz.")
    parser.add_argument("--order", type=int, default=4, help="IIR filter order.")
    parser.add_argument("--zerophase", type=int, default=1, choices=[0, 1], help="0 causal, 1 forward/backward.")
    parser.add_argument("--duration", type=float, default=12.0, help="Synthetic signal duration in seconds.")
    parser.add_argument("--seed", type=int, default=20260425, help="Random seed.")
    parser.add_argument("--ftype", default="butter", choices=["butter"], help="Python reference prototype.")
    parser.add_argument("--input-file", default="x.txt", help="Input file your FLIT executable reads.")
    parser.add_argument("--output-file", default="y.txt", help="Output file your FLIT executable writes.")
    parser.add_argument("--out", default="flit_vs_obspy_all_filters.png", help="Output figure.")
    parser.add_argument("--no-show", action="store_true", help="Save figures without opening GUI.")
    parser.add_argument(
        "--notch-mode",
        default="center_width",
        choices=["center_width", "bounds"],
        help=(
            "How to interpret notch flow/fhigh for Python reference. "
            "center_width: flow=center, fhigh=width. "
            "bounds: flow=low cutoff, fhigh=high cutoff."
        ),
    )

    args = parser.parse_args()

    zerophase = bool(args.zerophase)

    t, x = make_random_signal(args.fs, args.duration, args.seed)
    write_vec(args.input_file, x)

    cases = default_cases(args.fs)

    if args.notch_mode == "bounds":
        for c in cases:
            if c["type"] == "notch":
                center = c["flow"]
                width = c["fhigh"]
                c["flow"] = center - 0.5 * width
                c["fhigh"] = center + 0.5 * width
                c["label"] = f"notch/bandstop {c['flow']:.3g}-{c['fhigh']:.3g} Hz"

    results = []

    print()
    print("Validation summary")
    print("-" * 96)
    print(f"{'type':10s} {'fs':>8s} {'flow':>10s} {'fhigh':>10s} {'order':>6s} {'zero':>5s} {'max_abs':>16s} {'rel_l2':>16s}")
    print("-" * 96)

    for case in cases:
        y_fortran = run_fortran_exec(
            exe=args.exe,
            case=case,
            fs=args.fs,
            order=args.order,
            zerophase=zerophase,
            output_file=args.output_file,
        )

        if y_fortran.shape != x.shape:
            raise RuntimeError(
                f"Fortran output shape mismatch for {case['type']}: "
                f"expected {x.shape}, got {y_fortran.shape}."
            )

        ref_case = dict(case)
        if args.notch_mode == "bounds" and ref_case["type"] == "notch":
            ref_case["type"] = "bandstop"

        sos = scipy_reference_sos(
            case=ref_case,
            fs=args.fs,
            order=args.order,
            ftype=args.ftype,
        )
        y_ref = obspy_style_filter(x, sos, zerophase=zerophase)

        err = y_fortran - y_ref
        max_abs = float(np.max(np.abs(err)))
        rel_l2 = float(np.linalg.norm(err) / max(np.linalg.norm(y_ref), 1.0e-300))

        print(
            f"{case['type']:10s} "
            f"{args.fs:8.2f} "
            f"{case['flow']:10.3f} "
            f"{case['fhigh']:10.3f} "
            f"{args.order:6d} "
            f"{int(zerophase):5d} "
            f"{max_abs:16.6e} "
            f"{rel_l2:16.6e}"
        )

        results.append(
            {
                "case": case,
                "label": case["label"],
                "y_fortran": y_fortran,
                "y_ref": y_ref,
                "err": err,
                "max_abs": max_abs,
                "rel_l2": rel_l2,
            }
        )

    print("-" * 96)

    plot_all(
        t=t,
        x=x,
        results=results,
        fs=args.fs,
        out=args.out,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()