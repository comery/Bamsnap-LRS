import unittest
from src.bamsnap_lrs.reader import Read
from src.bamsnap_lrs.cigar import Segment
from src.bamsnap_lrs.pileup import base_pileup, pileup_to_pixels


class TestPileup(unittest.TestCase):
    def test_base_pileup(self):
        segs = [Segment("match", "=", 5, 5, 5)]
        r = Read(qname="q1", chrom="chr1", start=100, end=105, reverse=False, mapq=60, primary=True, supplementary=False, secondary=False, segments=segs, seq="AACGT")
        pile = base_pileup([r], 100, 105)
        self.assertEqual(pile[0]["A"], 1)
        self.assertEqual(pile[1]["A"], 1)
        self.assertEqual(pile[2]["C"], 1)
        self.assertEqual(pile[3]["G"], 1)
        self.assertEqual(pile[4]["T"], 1)
        bins = pileup_to_pixels(pile, 2)
        self.assertEqual(sum(v for k, v in bins[0].items() if k != "depth"), bins[0]["depth"])


if __name__ == "__main__":
    unittest.main()