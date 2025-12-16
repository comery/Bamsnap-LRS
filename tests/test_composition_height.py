import unittest
from src.bamsnap_lrs.reader import Read
from src.bamsnap_lrs.cigar import Segment
from src.bamsnap_lrs.renderer import render_snapshot


class TestCompositionHeight(unittest.TestCase):
    def test_depth_scaled_height(self):
        r1 = Read(qname="q1", chrom="chr1", start=0, end=1, reverse=False, mapq=60, primary=True, supplementary=False, secondary=False, segments=[Segment("match","=",1,1,1)], seq="A")
        reads = [r1]
        for i in range(10):
            reads.append(Read(qname=f"q2_{i}", chrom="chr1", start=1, end=2, reverse=False, mapq=60, primary=True, supplementary=False, secondary=False, segments=[Segment("match","=",1,1,1)], seq="A"))
        img = render_snapshot(reads, "chr1", 0, 2, width=2, read_height=6, detail="low", show_axis=False, show_composition=True, composition_height=20, style="default")
        px = img.load()
        h = img.size[1]
        bottom = 10 + 20
        def top_for_x(x):
            y = bottom
            while y > 0 and px[x, y-1] != (255,255,255):
                y -= 1
            return y
        t0 = top_for_x(0)
        t1 = top_for_x(1)
        self.assertLess(t1, t0)


if __name__ == "__main__":
    unittest.main()