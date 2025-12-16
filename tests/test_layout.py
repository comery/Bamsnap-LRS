import unittest
from src.bamsnap_lrs.cigar import Segment
from src.bamsnap_lrs.layout import segments_to_pixels


class TestLayout(unittest.TestCase):
    def test_segments_to_pixels(self):
        segs = [
            Segment("match", "=", 100, 100, 100),
            Segment("mismatch", "X", 1, 1, 1),
            Segment("ins", "I", 3, 0, 3),
            Segment("del", "D", 5, 5, 0),
            Segment("match", "=", 50, 50, 50),
        ]
        px = segments_to_pixels(segs, read_start=1000, region_start=1000, bp_per_px=10)
        self.assertTrue(any(t == "ins" for t, _, _ in px))
        lens = [x1 - x0 for t, x0, x1 in px if t != "ins"]
        self.assertTrue(sum(lens) > 0)


if __name__ == "__main__":
    unittest.main()