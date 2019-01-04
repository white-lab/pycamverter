from unittest import TestCase

from pycamv.fragment import ms_labels


class MSLabelsTest(TestCase):
    def test_labels(self):
        for key in ms_labels.LABEL_NUMBERS.keys():
            self.assertIn(key, ms_labels.LABEL_MZ_WINDOW)
            self.assertIn(key, ms_labels.LABEL_MASSES)
            self.assertIn(key, ms_labels.LABEL_NAMES)

        for key, val in ms_labels.LABEL_MASSES.items():
            self.assertLess(ms_labels.LABEL_MZ_WINDOW[key][0], min(val))

            # Skip greater check for TMT10plex
            if key in ["TMT10plex", "TMT11plex"]:
                continue

            self.assertGreater(ms_labels.LABEL_MZ_WINDOW[key][1], max(val))
