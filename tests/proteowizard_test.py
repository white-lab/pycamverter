from unittest import TestCase

from pycamv import proteowizard


class ProteoWizardTest(TestCase):
    def test_fetch_proteowizard(self):
        proteowizard.fetch_proteowizard()
