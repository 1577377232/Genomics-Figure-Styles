"""Synthetic-only checks for the compact profile; no project data."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
SPEC=importlib.util.spec_from_file_location("clean",Path(__file__).resolve().parents[1]/"python/manhattan_clean.py")
mod=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(mod)
class CompactProfileTest(unittest.TestCase):
    def test_inclusive_cutoff_unknown_and_full_point_counts(self):
        sv=pd.DataFrame({"chrom":["1","Z","1"],"position":[10,20,30],"P":[.01,1e-9,.5],"TE_coverage":[.5,.49,np.nan]})
        layout=pd.DataFrame({"chrom":["1","Z"],"length":[100,100]})
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/"panel.pdf"
            rec=mod.draw(np.array([1,40]),np.array([20,30]),np.array([.2,1e-10]),sv,layout,"Synthetic test",out)
            self.assertTrue(out.exists())
            self.assertEqual((rec["SNP_plotted"],rec["SV_plotted"],rec["TEV_plotted"]),(2,2,1))
            self.assertEqual(rec["SV_annotation_unknown"],1)
            self.assertEqual(rec["SNP_minP"],1e-10)
            self.assertGreater(rec["ymax"],10)
            self.assertIsNone(rec["threshold_P"])
            self.assertEqual(rec["outer_padding_bp"],12000000)
if __name__=="__main__": unittest.main()
