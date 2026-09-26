#!/usr/bin/env python3
"""
delly_analysis.py — Compare AA, SVRecalibrator, and Delly junction calls.

Reads consolidated.csv (one row per structural variant / breakpoint) and
produces homology/insertion concordance statistics plus a three-panel
length-distribution histogram (one panel each for AA, SVRecalibrator, Delly).

Validity conventions:
  * AA: ``homology_len`` empty = no result; 0 = valid blunt (exact ligation);
    negative = untemplated insertion (length = abs), positive = homology.
  * SVRecalibrator: same idea for ``sc_hom_len`` / ``sp_hom_len``. A valid
    SVRecalibrator result is a row where either column is populated.
  * Delly: ``delly_homlen`` missing = no result (0 = valid blunt);
    ``delly_inslen`` of 0 is invalid, same as missing (a real insertion call
    needs inslen > 0). Delly junction sequence is recovered from the split-read
    consensus: homology = consensus[consbp:consbp+homlen],
    insertion = consensus[consbp:consbp+inslen] (0-based ``delly_consbp``).

Only rows with ``orientation_match`` and ``delly_selected_by_svr_rows`` set
are treated as AA<->delly comparable (the delly call on other rows does not
correspond to the AA call). If ``ends_swapped`` is True, end 1 of AA
corresponds to end 2 of delly and vice versa.

Sequence matching is orientation-insensitive: a match is declared when the two
sequences are equal or reverse-complements of each other.

Usage:
    python analysis/svr_delly_analysis.py [path/to/consolidated.csv] [-o OUTPUT_PREFIX]

    The input defaults to consolidated.csv in the repository root, so the
    script works regardless of the current working directory. Outputs default
    to analysis/svr_delly_results/ (created if it does not exist).

Outputs:
    <prefix>_stats.txt   human-readable statistics (also printed to stdout)
    <prefix>_hist.png    stacked AA / SVRecalibrator / Delly length histograms
    <prefix>_hist.pdf    same figure as PDF
    <prefix>_part7.tsv   per-SV details for statistic 7 (delly valid, SVR invalid)
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def rev_comp(seq):
    """Reverse complement (tolerates lower case and N)."""
    if not isinstance(seq, str) or seq == "":
        return seq
    return "".join(COMP.get(b.upper(), "N") for b in reversed(seq))


def seqs_match(a, b):
    """Orientation-insensitive sequence equality (empty matches only empty)."""
    a = (a or "").upper()
    b = (b or "").upper()
    if a == "" or b == "":
        return a == b
    return a == b or rev_comp(a) == b


def parse_int_len(series):
    """Parse AA-style lengths ('2', '2.0', '-1.0') -> float Series with NaN."""
    return pd.to_numeric(series, errors="coerce")


def parse_strict_int(series):
    """Parse SVR/delly-style integer lengths -> float Series with NaN."""
    return pd.to_numeric(series, errors="coerce")


def delly_junction(df):
    """Return (signed_len, seq, cls) Series for delly rows.

    Validity (non-NaN signed_len) requires delly_homlen present. inslen > 0
    would indicate an insertion; inslen == 0 / missing means "no insertion".
    """
    homlen = parse_strict_int(df["delly_homlen"])
    inslen = parse_strict_int(df["delly_inslen"]).fillna(0)
    consbp = parse_strict_int(df["delly_consbp"])
    consensus = df["delly_consensus"].fillna("").astype(str)

    signed = pd.Series(np.nan, index=df.index, dtype=float)
    seq = pd.Series("", index=df.index, dtype=object)

    valid = homlen.notna()
    is_ins = valid & (inslen > 0)
    is_hom = valid & (~is_ins) & (homlen > 0)
    is_blunt = valid & (~is_ins) & (homlen == 0)

    signed[is_ins] = -inslen[is_ins]
    signed[is_hom] = homlen[is_hom]
    signed[is_blunt] = 0.0

    # Recover junction sequence from the consensus contig. Guard against
    # out-of-bounds slices (treat as unavailable -> empty -> non-match).
    for idx in df.index[is_ins | is_hom]:
        try:
            bp = int(consbp.loc[idx])
            ln = int(inslen.loc[idx]) if idx in df.index[is_ins] else int(homlen.loc[idx])
            con = consensus.loc[idx]
            if 0 <= bp and bp + ln <= len(con):
                seq.loc[idx] = con[bp : bp + ln]
        except (ValueError, TypeError):
            pass

    cls = pd.Series("invalid", index=df.index, dtype=object)
    cls[is_ins] = "insertion"
    cls[is_hom] = "homology"
    cls[is_blunt] = "blunt"
    return signed, seq, cls


def aa_junction(df):
    """Return (signed_len, seq, cls) Series for AA rows."""
    signed = parse_int_len(df["homology_len"])
    seq = df["homology_seq"].fillna("").astype(str)
    cls = pd.Series("invalid", index=df.index, dtype=object)
    cls[signed > 0] = "homology"
    cls[signed < 0] = "insertion"
    cls[signed == 0] = "blunt"
    return signed, seq, cls


def svr_signed_len(df):
    """Unified SVRecalibrator signed length.

    Uses sc when only sc is present, sp when only sp is present, and the
    (equal) value when both agree. Rows where sc and sp lengths disagree are
    NaN (excluded).
    """
    sc = parse_strict_int(df["sc_hom_len"])
    sp = parse_strict_int(df["sp_hom_len"])
    out = pd.Series(np.nan, index=df.index, dtype=float)
    only_sc = sc.notna() & sp.isna()
    only_sp = sp.notna() & sc.isna()
    both_agree = sc.notna() & sp.notna() & (sc == sp)
    out[only_sc] = sc[only_sc]
    out[only_sp] = sp[only_sp]
    out[both_agree] = sc[both_agree]
    return out, sc, sp


def positions_match_exactly(df):
    """Exact break-position agreement, honouring ends_swapped."""
    b1 = pd.to_numeric(df["break_pos1"], errors="coerce")
    b2 = pd.to_numeric(df["break_pos2"], errors="coerce")
    d1 = pd.to_numeric(df["delly_pos1"], errors="coerce")
    d2 = pd.to_numeric(df["delly_pos2"], errors="coerce")
    swapped = df["ends_swapped"].astype(str) == "True"
    straight = (~swapped) & b1.notna() & b2.notna() & d1.notna() & d2.notna()
    straight_match = straight & (b1 == d1) & (b2 == d2)
    swapped_match = swapped & b1.notna() & b2.notna() & d1.notna() & d2.notna()
    swapped_match = swapped_match & (b1 == d2) & (b2 == d1)
    return straight_match | swapped_match


def main(args_list=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_input = os.path.join(os.path.dirname(script_dir), "consolidated.csv")
    p.add_argument("input", nargs="?", default=default_input,
                   help="path to consolidated.csv "
                        "(default: consolidated.csv in the repository root)")
    default_prefix = os.path.join(script_dir, "svr_delly_results", "svr_delly_analysis")
    p.add_argument("-o", "--output-prefix", default=default_prefix,
                   help="prefix for stats/figure outputs "
                        "(default: analysis/svr_delly_results/svr_delly_analysis)")
    args = p.parse_args(args_list)

    out_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(args.input, dtype=str, keep_default_na=False)
    total = len(df)

    # ── comparable universe: delly call corresponds to the AA call ──
    comparable = (df["orientation_match"].astype(str) == "True") & (
        pd.to_numeric(df["delly_selected_by_svr_rows"], errors="coerce").fillna(0) != 0
    )

    aa_len, aa_seq, aa_cls = aa_junction(df)
    de_len, de_seq, de_cls = delly_junction(df)
    svr_len, sc_len, sp_len = svr_signed_len(df)

    aa_valid = aa_len.notna()
    de_valid = de_len.notna()
    svr_valid = sc_len.notna() | sp_len.notna()

    both_valid = comparable & aa_valid & de_valid
    n_both = int(both_valid.sum())

    def match_mask(cls):
        m = pd.Series(False, index=df.index)
        m[both_valid] = [
            (aa_cls[i] == cls)
            and (de_cls[i] == cls)
            and (abs(aa_len[i]) == abs(de_len[i]))
            and seqs_match(aa_seq[i], de_seq[i])
            for i in df.index[both_valid]
        ]
        return m

    m_overall = match_mask("homology") | match_mask("insertion") | match_mask("blunt")
    # blunt: no sequence to compare, length equality (0 == 0) suffices
    m_hom = match_mask("homology")
    m_ins = match_mask("insertion")
    m_blunt = both_valid & (aa_cls == "blunt") & (de_cls == "blunt")
    assert int(m_overall.sum()) == int(m_hom.sum() + m_ins.sum() + m_blunt.sum()), \
        "homology/insertion/blunt matches should partition overall matches"

    pos_ok = positions_match_exactly(df)
    m_overall_pos = m_overall & pos_ok
    m_hom_pos = m_hom & pos_ok
    m_ins_pos = m_ins & pos_ok
    m_blunt_pos = m_blunt & pos_ok

    n_delly_no_svr = int(((~svr_valid) & de_valid).sum())
    part7 = df.index[(~svr_valid) & de_valid]

    def compact(row, idx):
        def show(v):
            v = str(v)
            return v if v != "" else "-"

        junction = de_seq.loc[idx] if de_seq.loc[idx] != "" else "-"
        aa = (f"ends={row['break_chrom1']}:{row['break_pos1']}-"
              f"{row['break_chrom2']}:{row['break_pos2']} "
              f"svtype={show(row['break_sv_type'])} ori={row['break_orientation']} "
              f"hom_len={show(row['homology_len'])} seq={show(row['homology_seq'])}")
        de = (f"ends={show(row['delly_chrom1'])}:{show(row['delly_pos1'])}-"
              f"{show(row['delly_chrom2'])}:{show(row['delly_pos2'])} "
              f"svtype={show(row['delly_svtype'])} ori={show(row['delly_orientation'])} "
              f"homlen={show(row['delly_homlen'])} inslen={show(row['delly_inslen'])} "
              f"seq={junction}")
        return f"{row['svr_row_id']} | AA {aa} | delly {show(row['delly_id'])} {de}"

    # ── report ──
    def pct(n, d):
        return f"{100 * n / d:.1f}%" if d else "n/a"

    L = []
    L.append("Delly vs AA/SVRecalibrator junction comparison")
    L.append(f"Input: {args.input}  (n={total} SVs)")
    L.append(f"Comparable SVs (orientation_match + delly_selected_by_svr_rows): "
             f"{int(comparable.sum())}")
    L.append("")
    L.append("Validity:")
    L.append(f"  AA valid (homology_len present)               : {int(aa_valid.sum())}")
    L.append(f"  Delly valid (delly_homlen present)            : {int(de_valid.sum())}")
    L.append(f"  SVRecalibrator valid (sc or sp populated)     : {int(svr_valid.sum())}")
    L.append("")
    L.append("Statistics (comparable SVs only, unless noted):")
    L.append(f"  1. AA and delly both valid                    : {n_both} "
             f"({pct(n_both, int(comparable.sum()))} of comparable)")
    L.append(f"  2. AA/delly junction+length match (any class) : {int(m_overall.sum())} / {n_both} "
             f"({pct(int(m_overall.sum()), n_both)} of both-valid)")
    L.append(f"  3. AA/delly homology+length match             : {int(m_hom.sum())} / {n_both} "
             f"({pct(int(m_hom.sum()), n_both)})")
    L.append(f"  4. AA/delly insertion+length match            : {int(m_ins.sum())} / {n_both} "
             f"({pct(int(m_ins.sum()), n_both)})")
    L.append(f"  5. Both detect exact ligation (blunt)         : {int(m_blunt.sum())} / {n_both} "
             f"({pct(int(m_blunt.sum()), n_both)})")
    L.append("  6. Same as 2-5 but break positions match exactly (ends_swapped-aware):")
    L.append(f"     2'. overall match + exact positions        : {int(m_overall_pos.sum())}")
    L.append(f"     3'. homology match + exact positions       : {int(m_hom_pos.sum())}")
    L.append(f"     4'. insertion match + exact positions      : {int(m_ins_pos.sum())}")
    L.append(f"     5'. blunt match + exact positions          : {int(m_blunt_pos.sum())}")
    L.append(f"  7. Delly valid but SVRecalibrator invalid (all SVs): {n_delly_no_svr}")
    for idx in part7:
        L.append(f"       {compact(df.loc[idx], idx)}")
    L.append("")
    L.append("Histogram Inputs:")
    sc_sp_both = int((sc_len.notna() & sp_len.notna()).sum())
    sc_sp_disagree = int((sc_len.notna() & sp_len.notna() & (sc_len != sp_len)).sum())
    L.append(f"  SVR rows with both sc+sp results: {sc_sp_both}, "
             f"length disagreement (excluded): {sc_sp_disagree}")

    report = "\n".join(L)
    print(report)
    stats_path = args.output_prefix + "_stats.txt"
    with open(stats_path, "w") as f:
        f.write(report + "\n")
    print(f"Written: {stats_path}")

    # ── compact per-SV table for statistic 7 ──
    if len(part7):
        detail_cols = ["svr_row_id", "sample_id", "break_chrom1", "break_pos1",
                       "break_chrom2", "break_pos2", "break_orientation",
                       "homology_len", "homology_seq", "sc_pos1", "sc_pos2",
                       "sc_hom_len", "sc_hom", "sp_left_sv", "sp_right_sv",
                       "sp_hom_len", "hom", "delly_id", "delly_chrom1", "delly_pos1",
                       "delly_chrom2", "delly_pos2", "delly_orientation",
                       "delly_homlen", "delly_inslen", "delly_consbp",
                       "ends_swapped", "orientation_match",
                       "delly_selected_by_svr_rows"]
        details = df.loc[part7, [c for c in detail_cols if c in df.columns]].copy()
        details["delly_junction_seq"] = [de_seq.loc[i] for i in part7]
        details_path = args.output_prefix + "_part7.tsv"
        details.to_csv(details_path, sep="\t", index=False)
        print(f"Written: {details_path}")

    # ── histogram: one panel per caller, stacked vertically ──
    aa_vals = aa_len[aa_valid].astype(int)
    svr_vals = svr_len[svr_len.notna()].astype(int)
    # delly signed: -inslen if insertion else homlen (0 = blunt)
    de_vals = de_len[de_valid].astype(int)

    lo, hi = -50, 50
    bins = np.arange(lo - 0.5, hi + 1.5, 1.0)
    tracks = [
        (aa_vals, "AA", "#D65F5F"),
        (svr_vals, "SVRecalibrator", "#4878CF"),
        (de_vals, "Delly", "#6ACC65"),
    ]

    fig, axes = plt.subplots(len(tracks), 1, figsize=(9, 7), sharex=True)
    for ax, (vals, label, color) in zip(axes, tracks):
        n_clip_lo = int((vals < lo).sum())
        n_clip_hi = int((vals > hi).sum())
        ax.hist(vals.clip(lo, hi), bins=bins, color=color, alpha=0.78,
                edgecolor="white", linewidth=0.3)
        ax.axvline(0, color="black", linewidth=0.9, linestyle="--", zorder=3)
        ax.set_ylabel("SV count")
        ax.text(0.01, 0.93, f"{label} (n={len(vals)})", transform=ax.transAxes,
                ha="left", va="top", fontsize=10, fontweight="bold", color=color,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))
        notes = []
        if n_clip_lo:
            notes.append(f"+{n_clip_lo} in ≤{lo} bin")
        if n_clip_hi:
            notes.append(f"+{n_clip_hi} in ≥{hi} bin")
        if notes:
            ax.text(0.99, 0.93, "  ".join(notes), transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, color="gray")
    max_y = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, max_y)
    axes[-1].set_xlabel("Signed junction length (bp):  insertion (−)  |  blunt (0)  |  homology (+)")
    axes[0].set_title("Junction length distributions: AA vs SVRecalibrator vs Delly")
    # Label the ≤−50 and ≥50 endpoint bins explicitly so the cutoff is obvious.
    from matplotlib.ticker import FixedLocator
    bot = axes[-1]
    ticks = sorted(set(list(bot.get_xticks()) + [lo, hi]))
    bot.xaxis.set_major_locator(FixedLocator(ticks))
    bot.set_xticklabels(
        ["≤−50" if t == lo else "≥50" if t == hi else str(int(t)) for t in ticks],
    )
    fig.tight_layout()

    for ext in (".png", ".pdf"):
        path = args.output_prefix + "_hist" + ext
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Written: {path}")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
