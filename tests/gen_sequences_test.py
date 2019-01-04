from unittest import TestCase

from pycamv.fragment import gen_sequences


class GenSequencesTest(TestCase):
    def test_gen_possible_seq_test(self):
        seqs = list(
            gen_sequences.gen_possible_seq(
                "IEFTTER",
                [(1, "Phospho", ["T"])],
            )
        )
        self.assertEqual(len(seqs), 2)
        self.assertIn(
            [
                ("N-term", ()),
                ("I", ()),
                ("E", ()),
                ("F", ()),
                ("T", ("Phospho",)),
                ("T", ()),
                ("E", ()),
                ("R", ()),
                ("C-term", ()),
            ],
            seqs,
        )
        self.assertIn(
            [
                ("N-term", ()),
                ("I", ()),
                ("E", ()),
                ("F", ()),
                ("T", ()),
                ("T", ("Phospho",)),
                ("E", ()),
                ("R", ()),
                ("C-term", ()),
            ],
            seqs,
        )
