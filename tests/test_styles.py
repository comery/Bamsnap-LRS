import unittest
from src.bamsnap_lrs.styles import color_for_base, color_for_type, shade_by_mapq


class TestStyles(unittest.TestCase):
    def test_base_colors(self):
        self.assertEqual(color_for_base('A'), (80,160,80))
        self.assertEqual(color_for_base('C'), (80,120,200))
        self.assertEqual(color_for_base('G'), (255,165,0))
        self.assertEqual(color_for_base('T'), (220,90,90))
        self.assertEqual(color_for_base('N'), (160,160,160))

    def test_type_colors(self):
        self.assertEqual(color_for_type('match'), (200,200,200))
        self.assertEqual(color_for_type('ins'), (60,120,200))

    def test_mapq_shade(self):
        c = (200,200,200)
        self.assertEqual(shade_by_mapq(c, 0), (40,40,40))
        self.assertEqual(shade_by_mapq(c, 60), (200,200,200))


if __name__ == '__main__':
    unittest.main()