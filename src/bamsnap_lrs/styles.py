from typing import Tuple


BASE_COLORS = {
    "A": (80, 160, 80),
    "C": (80, 120, 200),
    "G": (255, 165, 0),
    "T": (220, 90, 90),
    "N": (160, 160, 160),
}


TYPE_COLORS = {
    "match": (180, 180, 180),  # Gray for matches (JBrowse style)
    "mismatch": (200, 60, 60),
    "ins": (128, 0, 128),  # Purple for insertions (like JBrowse)
    "del": (128, 128, 128),  # Dark gray for deletions
    "ref_skip": (240, 240, 240),
    "soft": (171, 219, 227),  # Gray like match
    "hard": (118, 181, 197),  # Gray like match
}

# JBrowse-style colors for mismatches (more vibrant)
# Used for read mismatch and variant bases in coverage track
MISMATCH_COLORS = {
    "A": (255, 0, 0),      # Red
    "C": (0, 0, 255),      # Blue
    "G": (0, 128, 0),      # Green
    "T": (255, 165, 0),    # Orange
    "N": (128, 128, 128),  # Gray
}


STRAND_COLORS = {
    "+": (170, 170, 170),
    "-": (120, 120, 160),
}


def shade_by_mapq(color: Tuple[int, int, int], mapq: int) -> Tuple[int, int, int]:
    f = min(1.0, max(0.2, mapq / 60.0))
    return tuple(int(c * f) for c in color)


def color_for_base(base: str) -> Tuple[int, int, int]:
    return BASE_COLORS.get(base.upper(), BASE_COLORS["N"])


def color_for_type(t: str) -> Tuple[int, int, int]:
    return TYPE_COLORS.get(t, (100, 100, 100))
