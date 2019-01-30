from queue import Queue
import tempfile
from unittest import TestCase

from pycamv import export, fragment, scan, search


class ExportTest(TestCase):
    def test_export_to_sql(self):
        fd, path = tempfile.mkstemp(suffix='.db')
        filename = "dummy/path.msf"

        pep_query = search.pep_query.PeptideQuery(
            accessions=('Fgr',),
            prot_descs=('Tyrosine-protein kinase Fgr',),
            uniprot_accessions=('P09769',),
            full_seqs=(
                "MGCVFCKKLEPVATAKEDAGLEGDFRSYGAADHYGPDPTKARPASSFAHIPNYSNFSSQA"
                "INPGFLDSGTIRGVSGIGVTLFIALYDYEARTEDDLTFTKGEKFHILNNTEGDWWEARSL"
                "SSGKTGCIPSNYVAPVDSIQAEEWYFGKIGRKDAERQLLSPGNPQGAFLIRESETTKGAY"
                "SLSIRDWDQTRGDHVKHYKIRKLDMGGYYITTRVQFNSVQELVQHYMEVNDGLCNLLIAP"
                "CTIMKPQTLGLAKDAWEISRSSITLERRLGTGCFGDVWLGTWNGSTKVAVKTLKPGTMSP"
                "KAFLEEAQVMKLLRHDKLVQLYAVVSEEPIYIVTEFMCHGSLLDFLKNPEGQDLRLPQLV"
                "DMAAQVAEGMAYMERMNYIHRDLRAANILVGERLACKIADFGLARLIKDDEYNPCQGSKF"
                "PIKWTAPEAALFGRFTIKSDVWSFGILLTELITKGRIPYPGMNKREVLEQVEQGYHMPCP"
                "PGCPASLYEAMEQTWRLDPEERPTFEYLQSFLEDYFTSAEPQYQPGDQT",
            ),
            query=1,
            filename=filename,
            pep_score=50,
            pep_exp_mz=611.23,
            pep_exp_z=2,
            pep_seq="GAYSLSIR",
            pep_var_mods=(
                (1, "Phospho", ("Y",)),
            ),
            pep_fixed_mods=(
                (1, "TMT6plex", ("N-term",)),
            ),
            scan=1234,
            rank_pos=1,
            quant_scan=1235,
        )
        scan_query = scan.scans.ScanQuery(
            scan,
            precursor_scan=1221,
            window_offset=[.05, .05],
            isolation_mz=611.231,
            collision_type='HCD',
            c13_num=1,
            basename=pep_query.basename,
        )
        mh_peak = fragment.compare.PeptideHit(
            611.231,
            1e5,
            name='MH^{+}',
            score=1.2,
            predicted_mz=611.23,
            match_list=[],
        )
        peaks = [
            mh_peak,
        ]
        queue = Queue()
        queue.put((
            pep_query,
            (
                ("N-term", ("TMT6plex",)),
                ("G", ()),
                ("A", ()),
                ("Y", ("Phospho",)),
                ("S", ()),
                ("L", ()),
                ("S", ()),
                ("I", ()),
                ("R", ()),
                ("C-term", ()),
            ),
            "maybe",
            peaks,
            [mh_peak],
            [mh_peak],
        ))
        scan_mapping = {pep_query: scan_query}

        export.export.export_to_sql(
            path,
            queue,
            scan_mapping,
            filename,
            [],
            total_num_seq=1,
        )
