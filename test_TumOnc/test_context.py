import unittest
import TumOnc.context as ctx
import pandas as pd

class test_context(unittest.TestCase):
    def test_ctx3id(self):
        for c in ctx.ctx3_id:
            self.assertEqual(ctx.id_ctx3[ctx.ctx3_id[c]],c)
        self.assertEqual(ctx.ctx3_id['ACG'],3+4*2+16*0)

    def test_ctx5id(self):
        for c in ctx.ctx5_id:
            self.assertEqual(ctx.id_ctx5[ctx.ctx5_id[c]],c)
        self.assertEqual(ctx.ctx5_id['TACGT'],(3+4*2+16*0)*16+1+1*4)

    def test_get_ctx_bins(self):
        td = pd.DataFrame(['AGCTT','AGCTT','AGCTT','AGCTT',
                           'AGTGA','AGTGA',
                           'GCCGG',
                           'ACCGA','ACCGT','ACCGG','ACCGC'],columns=['ctx5'])
        td['ctx3'] = td['ctx5'].str[1:4]
        bins = ctx.get_ctx_bins(td)
        # print(td)
        # print(bins)
        self.assertEqual(bins.loc['AGCTT','c5id'],0)
        self.assertEqual(bins.loc['AGTGA','c5id'],1)
        self.assertGreaterEqual(bins.loc['GCCGG','c5id'],2)
        self.assertLessEqual(bins.loc['GCCGG','c5id'],6)
        self.assertGreater(bins.loc['GCCGT','c5id'],6)
        for c1 in 'AGCT':
            for c2 in 'AGCT':
                self.assertEqual(bins.loc[c1+'CCG'+c2,'c3id'],0)
                self.assertEqual(bins.loc[c1+'GCT'+c2,'c3id'],1)
            
if __name__ == '__main__':
    unittest.main()
