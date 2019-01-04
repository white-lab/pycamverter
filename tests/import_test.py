from unittest import TestCase


class ImportTest(TestCase):
    def test_imports(self):
        from pycamv import (
            main, proteowizard, regexes, utils, version,
        )
        from pycamv import camv_mat
        from pycamv.export import (
            export, migrations, sql,
        )
        from pycamv.fragment import (
            compare, fragments, gen_sequences, losses, masses, ms_labels,
            validate,
        )
        from pycamv.scan import (
            scan_list, scans,
        )
        from pycamv.search import (
            discoverer, mascot, pep_query, search,
        )
