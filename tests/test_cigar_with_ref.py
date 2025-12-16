import unittest
from src.bamsnap_lrs.cigar import from_cigar_with_ref


class TestCigarWithRef(unittest.TestCase):
    def test_refine_m(self):
        read = "ACGTACGT"
        ref = "ACGTTCGT"
        segs = from_cigar_with_ref("8M", read, ref)
        types = [(s.type, s.length) for s in segs]
        self.assertEqual(types, [("match", 4), ("mismatch", 1), ("match", 3)])


if __name__ == "__main__":
    unittest.main()