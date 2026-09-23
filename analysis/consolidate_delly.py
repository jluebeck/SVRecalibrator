#!/usr/bin/env python3
"""Layer single-sample Delly VCF calls onto archived SVRecalibrator tables."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

CT_ORIENTATION = {'3to5': '+-', '5to3': '-+', '3to3': '++', '5to5': '--'}


def read_table(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames or [], list(reader)


def write_csv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_vcf(path):
    """Preserve raw INFO/FORMAT as well as fields used for matching. No filtering."""
    calls = []
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#CHROM'):
                if len(fields) != 10:
                    raise ValueError(f'{path}: expected a single-sample VCF')
                vcf_sample = fields[9]
                continue
            if line.startswith('#') or not line.strip():
                continue
            info = dict(item.split('=', 1) if '=' in item else (item, True)
                        for item in fields[7].split(';'))
            kind = info.get('SVTYPE', '')
            end = info.get('POS2') if kind in ('BND', 'TRA') else info.get('END')
            if end is None:
                raise ValueError(f'{path}, {fields[2]}: missing partner coordinate')
            call = dict(delly_key=f'{path.name}:{len(calls)+1}', delly_id=fields[2],
                        delly_vcf_sample=vcf_sample, delly_chrom1=fields[0],
                        delly_pos1=int(fields[1]), delly_chrom2=info.get('CHR2', fields[0]),
                        delly_pos2=int(end), delly_orientation=CT_ORIENTATION.get(info.get('CT'), ''),
                        delly_filter=fields[6], delly_qual=fields[5], delly_alt=fields[4],
                        delly_precise='PRECISE' in info, delly_info=fields[7],
                        delly_format=fields[8], delly_genotype=fields[9])
            for key in ['SVTYPE', 'CT', 'HOMLEN', 'INSLEN', 'CONSENSUS', 'CONSBP',
                        'PE', 'SR', 'MAPQ', 'SRMAPQ', 'CIPOS', 'CIEND', 'SVMETHOD']:
                call['delly_' + key.lower()] = info.get(key, '')
            calls.append(call)
    return calls


def chrom(value):
    return value.removeprefix('chr')


def candidates(row, calls, window):
    """Match both original AA ends; swapping ends also swaps their signs."""
    found = []
    for call in calls:
        alternatives = []
        for swapped in (False, True):
            a, b = (2, 1) if swapped else (1, 2)
            if (chrom(row['break_chrom1']), chrom(row['break_chrom2'])) != (
                    chrom(call[f'delly_chrom{a}']), chrom(call[f'delly_chrom{b}'])):
                continue
            d1 = call[f'delly_pos{a}'] - int(row['break_pos1'])
            d2 = call[f'delly_pos{b}'] - int(row['break_pos2'])
            if max(abs(d1), abs(d2)) > window:
                continue
            ori = call['delly_orientation'][::-1] if swapped else call['delly_orientation']
            match = bool(ori) and ori == row['break_orientation']
            alternatives.append((not match, max(abs(d1), abs(d2)), abs(d1)+abs(d2),
                                 swapped, d1, d2))
        if alternatives:
            mismatch, distance, total, swapped, d1, d2 = min(alternatives)
            found.append(dict(call, ends_swapped=swapped, orientation_match=not mismatch,
                              delta_pos1=d1, delta_pos2=d2, max_distance=distance,
                              total_distance=total))
    return sorted(found, key=lambda c: (not c['orientation_match'], c['max_distance'],
                                        c['total_distance'], c['delly_key']))


def consolidate(sv_dir, delly_dir, output, window):
    if not sv_dir.is_dir() or not delly_dir.is_dir():
        raise ValueError('Both input directories must exist')
    tables = sorted(sv_dir.glob('*/final_augmented.tsv'))
    if not tables:
        raise ValueError('No final_augmented.tsv files found')
    output.mkdir(parents=True, exist_ok=True)
    all_rows, all_candidates, inventory, all_delly, inputs = [], [], [], [], []
    used = defaultdict(int)
    for directory in sorted(p for p in sv_dir.iterdir() if p.is_dir()):
        path = directory / 'final_augmented.tsv'
        sample = directory.name
        vcf = delly_dir / (sample + '.vcf')
        if not vcf.exists():
            vcf = delly_dir / (sample + '.vcf.gz')
        calls = read_vcf(vcf) if vcf.exists() else []
        fields, rows = read_table(path) if path.exists() else ([], [])
        if rows and not {'break_chrom1', 'break_pos1', 'break_chrom2', 'break_pos2', 'break_orientation'} <= set(fields):
            raise ValueError(f'{path}: required breakpoint columns missing')
        inventory.append(dict(sample_id=sample, svr_status='missing' if not path.exists() else 'nonempty' if rows else 'empty',
                              svr_rows=len(rows), done_flag=(directory/'done.flag').exists(),
                              delly_status='present' if vcf.exists() else 'missing', delly_calls=len(calls)))
        for source in [path, vcf]:
            if source.exists():
                inputs.append(dict(path=str(source.resolve()), sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
        for index, row in enumerate(rows, 1):
            identity = dict(sample_id=sample, svr_row_id=f'{sample}:{index}', source_row=index,
                            source_table=str(path.resolve()))
            nearby = candidates(row, calls, window)
            compatible = [c for c in nearby if c['orientation_match']]
            result = dict(identity, **row, candidate_count=len(nearby), compatible_count=len(compatible),
                          match_status='matched' if compatible else 'orientation_discordant' if nearby else
                          'no_candidate' if vcf.exists() else 'missing_vcf')
            if compatible:
                best = compatible[0]
                result.update(best)
                result['equally_near_count'] = sum((c['max_distance'], c['total_distance']) ==
                                                  (best['max_distance'], best['total_distance']) for c in compatible)
                used[(sample, best['delly_key'])] += 1
            all_rows.append(result)
            all_candidates.extend(dict(identity, **c) for c in nearby)
        all_delly.extend(dict(sample_id=sample, **c) for c in calls)
    for row in all_rows:
        row['delly_selected_by_svr_rows'] = used[(row['sample_id'], row.get('delly_key', ''))]
    for call in all_delly:
        call['selected_by_svr_rows'] = used[(call['sample_id'], call['delly_key'])]
    write_csv(output/'consolidated.csv', all_rows, None if all_rows else ['sample_id', 'svr_row_id', 'match_status'])
    write_csv(output/'candidates.csv', all_candidates, None if all_candidates else ['sample_id', 'svr_row_id', 'delly_key'])
    write_csv(output/'sample_inventory.csv', inventory)
    write_csv(output/'delly_calls.csv', all_delly, None if all_delly else ['sample_id', 'delly_key'])
    represented = {r['sample_id'] for r in inventory}
    extra = sorted(p.name for p in delly_dir.glob('*.vcf*') if p.name.removesuffix('.gz').removesuffix('.vcf') not in represented)
    (output/'manifest.json').write_text(json.dumps(dict(window_bp=window, anchor='original break_pos1/2; 1-based as stored',
        matching='both ends within window, CT orientation compatible; nearest max then summed distance; ties by record key',
        filters='all records retained, including LowQual and reference genotypes',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        delly_only_files=extra, inputs=inputs), indent=2)+'\n')
    print(f'Wrote {len(all_rows)} SV rows and {len(all_candidates)} candidate pairs to {output}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sv_dir', type=Path)
    parser.add_argument('delly_dir', type=Path)
    parser.add_argument('-o', '--outdir', type=Path, required=True)
    parser.add_argument('--window', type=int, default=100, help='Maximum distance at each original breakpoint (default 100 bp)')
    args = parser.parse_args()
    if args.window < 0:
        parser.error('--window must be nonnegative')
    consolidate(args.sv_dir, args.delly_dir, args.outdir, args.window)


if __name__ == '__main__':
    main()
