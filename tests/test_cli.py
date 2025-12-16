import io
import json
import sys
import unittest

from src.bamsnap_lrs.cli import main


class TestCLI(unittest.TestCase):
    def test_parse_cigar_cli(self):
        argv = ["prog", "parse-cigar", "--cigar", "5=2X3I4D"]
        backup = sys.argv
        sys.argv = argv
        buf = io.StringIO()
        backup_stdout = sys.stdout
        sys.stdout = buf
        try:
            main()
        finally:
            sys.argv = backup
            sys.stdout = backup_stdout
        out = buf.getvalue().strip()
        data = json.loads(out)
        self.assertTrue(len(data) >= 3)


if __name__ == "__main__":
    unittest.main()