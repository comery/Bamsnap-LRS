import unittest
from src.bamsnap_lrs.cigar import parse_cigar_string, from_cigar_md_cs, Segment


class TestCigar(unittest.TestCase):
    def test_parse_cigar_string(self):
        self.assertEqual(parse_cigar_string("10M"), [("M", 10)])
        self.assertEqual(parse_cigar_string("5=5X"), [("=", 5), ("X", 5)])

    def test_basic_ops(self):
        segs = from_cigar_md_cs("3S10M2S")
        self.assertEqual([(s.type, s.length) for s in segs], [("soft", 3), ("match", 10), ("soft", 2)])

    def test_md_split(self):
        segs = from_cigar_md_cs("10M", md="5A5")
        self.assertEqual([(s.type, s.length) for s in segs], [("match", 5), ("mismatch", 1), ("match", 4)])

    def test_indels(self):
        segs = from_cigar_md_cs("5M2I3M1D2M", md="5A3^C2")
        types = [s.type for s in segs]
        self.assertIn("ins", types)
        self.assertIn("del", types)

    def test_cs(self):
        segs = from_cigar_md_cs("10M", cs=":5*at+gg:5")
        self.assertEqual([(s.type, s.length) for s in segs], [("match", 5), ("mismatch", 1), ("ins", 2), ("match", 5)])


if __name__ == "__main__":
    unittest.main()