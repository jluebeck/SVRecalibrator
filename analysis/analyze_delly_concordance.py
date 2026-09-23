#!/usr/bin/env python3
"""Informal, row-level concordance summary, kept separate from consolidation."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from consolidate_delly import write_csv


def number(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


def sequence_relation(a, b):
    a, b = a.upper(), b.upper()
    if not a or not b or set(a+b) - set('ACGT'):
        return 'unavailable'
    reverse = b.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    if a == b:
        return 'exact'
    return 'reverse_complement' if a == reverse else 'different'


def analyze(path, output):
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    matched = [r for r in rows if r['match_status'] == 'matched']
    summary = dict(svr_rows=len(rows), match_status=dict(Counter(r['match_status'] for r in rows)),
                   matched_filters=dict(Counter(r['delly_filter'] for r in matched)),
                   ambiguous_nearest_rows=sum(number(r.get('equally_near_count')) > 1 for r in matched),
                   reused_delly_rows=sum(number(r['delly_selected_by_svr_rows']) > 1 for r in matched))
    details = []
    for row in matched:
        detail = {key: row[key] for key in ['sample_id', 'svr_row_id', 'delly_key', 'delly_filter', 'delly_precise']}
        for method, left, right, hom in [('split', 'sp_left_sv', 'sp_right_sv', 'sp_hom_len'),
                                          ('scaffold', 'sc_pos1', 'sc_pos2', 'sc_hom_len')]:
            swapped = row['ends_swapped'] == 'True'
            a, b = (2, 1) if swapped else (1, 2)
            p1, p2 = number(row.get(left)), number(row.get(right))
            distance = None if p1 is None or p2 is None else max(
                abs(p1-int(row[f'delly_pos{a}'])), abs(p2-int(row[f'delly_pos{b}'])))
            detail[method+'_max_distance'] = distance
            hl, dl = number(row.get(hom)), number(row.get('delly_homlen'))
            # Negative SVR values are insertions; absent Delly HOMLEN is unknown.
            detail[method+'_homlen_relation'] = ('unavailable' if hl is None or dl is None else
                'svr_insertion' if hl < 0 else 'equal' if hl == dl else 'different')
        detail['split_scaffold_sequence_relation'] = sequence_relation(row.get('hom', ''), row.get('sc_hom', ''))
        details.append(detail)
    for label, subset in [('all_matched', details), ('pass_precise', [r for r in details if r['delly_filter'] == 'PASS' and r['delly_precise'] == 'True'])]:
        stats = {'rows': len(subset)}
        for method in ['split', 'scaffold']:
            distances = [r[method+'_max_distance'] for r in subset if r[method+'_max_distance'] is not None]
            stats[method] = dict(comparable_coordinates=len(distances), within_bp={str(t):sum(d <= t for d in distances) for t in [0, 1, 5, 10, 50, 100]},
                                 homlen_relation=dict(Counter(r[method+'_homlen_relation'] for r in subset)))
        stats['split_scaffold_sequence_relation'] = dict(Counter(r['split_scaffold_sequence_relation'] for r in subset))
        summary[label] = stats
    summary['limitations'] = [
        'Informal SVRecalibrator-ascertained row-level comparison; not precision/recall or a truth set.',
        'Matching is anchored to original AA coordinates; nearby orientation-compatible calls need not be identical events.',
        'Raw reported coordinates are compared without base-offset or microhomology-shift normalization.',
        'Delly HOMLEN permits edits and is not necessarily equivalent to exact SVRecalibrator homology.',
        'Delly CONSENSUS is a full junction consensus, not a homology sequence. No direct homology-sequence comparison is inferred.',
        'Reverse complements here compare split/scaffold sequences only, within matched rows; short motifs alone do not prove equivalent junctions.',
        'Blank and N/A lengths remain unknown; negative SVRecalibrator lengths mean insertion.',
    ]
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output/'concordance_details.csv', details, None if details else ['sample_id', 'svr_row_id'])
    (output/'concordance_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('consolidated', type=Path)
    parser.add_argument('-o', '--outdir', type=Path, required=True)
    args = parser.parse_args()
    analyze(args.consolidated, args.outdir)


if __name__ == '__main__':
    main()
