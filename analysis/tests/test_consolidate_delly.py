import tempfile
import unittest
from pathlib import Path
from consolidate_delly import candidates, read_vcf, consolidate
from analyze_delly_concordance import sequence_relation, number


class DellyTests(unittest.TestCase):
    def test_swapped_ends_keep_breakend_signs(self):
        row = dict(break_chrom1='chr1', break_pos1='100', break_chrom2='chr2',
                   break_pos2='200', break_orientation='+-')
        call = dict(delly_chrom1='2', delly_pos1=201, delly_chrom2='1',
                    delly_pos2=99, delly_orientation='-+', delly_key='x')
        hit = candidates(row, [call], 1)[0]
        self.assertTrue(hit['ends_swapped'])
        self.assertTrue(hit['orientation_match'])
        self.assertEqual((hit['delta_pos1'], hit['delta_pos2']), (-1, 1))
        self.assertEqual(candidates(row, [call], 0), [])
        call['delly_orientation'] = '++'
        self.assertFalse(candidates(row, [call], 1)[0]['orientation_match'])

    def test_bnd_uses_pos2_not_end_and_preserves_lowqual(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'sample.vcf'
            path.write_text('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
                            'chr1\t100\tbnd\tA\tA[chr2:200[\t10\tLowQual\tSVTYPE=BND;END=101;CHR2=chr2;POS2=200;CT=3to5;PRECISE;HOMLEN=0\tGT\t0/1\n')
            call = read_vcf(path)[0]
            self.assertEqual(call['delly_pos2'], 200)
            self.assertEqual(call['delly_orientation'], '+-')
            self.assertEqual(call['delly_homlen'], '0')
            self.assertEqual(call['delly_filter'], 'LowQual')
            self.assertEqual(call['delly_inslen'], '')

    def test_sequences_and_missing_lengths(self):
        self.assertEqual(sequence_relation('AAC', 'GTT'), 'reverse_complement')
        self.assertEqual(sequence_relation('AT', 'AT'), 'exact')
        self.assertEqual(sequence_relation('', ''), 'unavailable')
        self.assertEqual(sequence_relation('N', 'N'), 'unavailable')
        self.assertIsNone(number('N/A'))
        self.assertEqual(number('-3'), -3)

    def test_compatible_candidate_beats_closer_orientation_conflict(self):
        row = dict(break_chrom1='chr1', break_pos1='100', break_chrom2='chr1',
                   break_pos2='500', break_orientation='+-')
        close = dict(delly_chrom1='chr1', delly_pos1=100, delly_chrom2='chr1',
                     delly_pos2=500, delly_orientation='++', delly_key='close')
        farther = dict(close, delly_pos1=105, delly_pos2=505,
                       delly_orientation='+-', delly_key='farther')
        hits = candidates(row, [close, farther], 5)
        self.assertEqual([h['delly_key'] for h in hits], ['farther', 'close'])
        self.assertEqual(len(candidates(row, [farther], 4)), 0)

    def test_empty_and_missing_outputs_are_inventory_entries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); sv = root/'sv'; vcf = root/'vcf'; out = root/'out'
            (sv/'empty').mkdir(parents=True); (sv/'missing').mkdir(); vcf.mkdir()
            (sv/'empty'/'final_augmented.tsv').write_text('\n')
            consolidate(sv, vcf, out, 100)
            import csv
            with (out/'sample_inventory.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([r['svr_status'] for r in rows], ['empty', 'missing'])
            self.assertTrue(all(r['delly_status'] == 'missing' for r in rows))


if __name__ == '__main__':
    unittest.main()
