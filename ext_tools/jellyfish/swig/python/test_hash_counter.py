#  This file is part of Jellyfish.
#
#  This work is dual-licensed under 3-Clause BSD License or GPL 3.0.
#  You can choose between one of them if you use this work.
#
#  `SPDX-License-Identifier: BSD-3-Clause OR  GPL-3.0`

import unittest
import sys
import random


import dna_jellyfish as jf

class TestHashCounter(unittest.TestCase):
    def setUp(self):
        jf.MerDNA.k(100)
        self.hash = jf.HashCounter(1024, 5)

    def test_info(self):
        self.assertEqual(100, jf.MerDNA.k())
        self.assertEqual(1024, self.hash.size())
        self.assertEqual(5, self.hash.val_len())

    def test_add(self):
        mer  = jf.MerDNA()
        good = True
        for i in range(1000):
            mer.randomize()
            val = random.randrange(1000)
            good = good and self.hash.add(mer, val)
            if not good: break
            if i % 3 > 0:
                nval = random.randrange(1000)
                val  = val + nval
                if i % 3 == 1:
                    good = good and (not self.hash.add(mer, nval))
                else:
                    good = good and self.hash.update_add(mer, nval)
            if not good: break
            good = good and (val == self.hash.get(mer)) and (val == self.hash[mer])
            if not good: break
        self.assertTrue(good)


if __name__ == '__main__':
    data = sys.argv.pop(1)
    unittest.main()
