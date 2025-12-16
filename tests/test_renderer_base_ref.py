import unittest
from src.bamsnap_lrs.reader import Read
from src.bamsnap_lrs.cigar import Segment
from src.bamsnap_lrs.renderer import render_snapshot
from src.bamsnap_lrs.styles import color_for_base, color_for_type


class TestRendererBaseRef(unittest.TestCase):
    def test_base_gray_when_match_ref(self):
        segs = [Segment("match", "=", 4, 4, 4)]
        r = Read(qname="q1", chrom="chr1", start=0, end=4, reverse=False, mapq=60, primary=True, supplementary=False, secondary=False, segments=segs, seq="ACGT")
        img = render_snapshot([r], "chr1", 0, 4, width=4, read_height=4, detail="high", show_axis=False, show_composition=False, color_by="base", ref_seq="ACGT", style="default")
        px = img.load()
        gray = color_for_type("match")
        self.assertEqual(px[0, 1], gray)
        self.assertEqual(px[1, 1], gray)
        self.assertEqual(px[2, 1], gray)
        self.assertEqual(px[3, 1], gray)

    def test_base_colored_when_mismatch_ref(self):
        segs = [Segment("match", "=", 4, 4, 4)]
        r = Read(qname="q1", chrom="chr1", start=0, end=4, reverse=False, mapq=60, primary=True, supplementary=False, secondary=False, segments=segs, seq="ACCT")
        img = render_snapshot([r], "chr1", 0, 4, width=4, read_height=4, detail="high", show_axis=False, show_composition=False, color_by="base", ref_seq="ACGT", style="default")
        px = img.load()
        gray = color_for_type("match")
        self.assertEqual(px[0, 1], gray)
        self.assertEqual(px[1, 1], gray)
        self.assertEqual(px[2, 1], color_for_base("C"))
        self.assertEqual(px[3, 1], gray)


if __name__ == "__main__":
    unittest.main()