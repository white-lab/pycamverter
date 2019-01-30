import hashlib
import os
import requests
import tempfile
from unittest import TestCase

try:
    from tempfile import TemporaryDirectory
except ImportError:
    from pycamv.utils import TemporaryDirectory

from pycamv import search, fragment

PD_URL_BASE = (
    "https://media.githubusercontent.com/media/"
    "white-lab/pycamverter-data/master/"
)

PD14_URL_BASE = 'Test%20PD1.4/'
PD22_URL_BASE = 'Test%20PD2.2/'

PD14_URLS = [
    ("CK-C1-pY.msf", "96eda5b0e47f615cf000d1c5d3ecc8cd"),
    ("FAD-H1-Global.msf", "497ba1841faac883619d4f91a86a95cc"),
    ("FAD-H1-MPM2.msf", "5def08356bcfa5b679835e4a23dd1396"),
    ("Tau-4moHL1-Global.msf", "95d8089b5e4657b348bea0868c655478"),
]
PD14_RAWS = {
    "CK-C1-pY.msf": [
        (
            "2015-11-13-CKC1-pY-imac30-elute-pre50-col40.raw",
            "e10ba1aeb6ae587549aa40e327eace93",
        ),
    ],
    "FAD-H1-Global.msf": [
        (
            "2017-02-27-FADH1-pY-sup10-pre48-col99.RAW",
            "e4b02c27fb5030a60ecaced28c56db66",
        ),
    ],
    "FAD-H1-MPM2.msf": [
        (
            "2017-03-16-FADH1-MPM2-NTA-pre100-col96.raw",
            "ca28737c529cbe871fa4095b4c2496f1",
        ),
    ],
    "Tau-4moHL1-Global.msf": [
        (
            "2017-10-13-TauHL1-pY-sup10-pre126-col113.RAW",
            "e30d3bcbbbe1cd83ca387e721039a4c0",
        ),
    ],
}

PD22_URLS = [
]


class ValidateTest(TestCase):
    def fetch_url(self, url, dir, md5hash=None):
        response = requests.get(url, stream=True)
        path = os.path.join(dir, url.split('/')[-1])
        hash_md5 = hashlib.md5()

        with open(path, 'wb') as f:
            for block in response.iter_content(1024):
                hash_md5.update(block)
                f.write(block)

        if md5hash is not None and hash_md5.hexdigest() != md5hash:
            raise Exception(
                "MD5 hash of {} does not match record: {} != {}"
                .format(url, md5hash, hash_md5.hexdigest())
            )

        return path

    def test_validate_pd14(self):
        for url, md5 in PD14_URLS:
            with TemporaryDirectory() as tmpdir:
                search_path = self.fetch_url(
                    PD_URL_BASE + PD14_URL_BASE + url,
                    tmpdir,
                    md5hash=md5,
                )
                raw_paths = [
                    self.fetch_url(
                        PD_URL_BASE + PD14_URL_BASE + raw,
                        tmpdir,
                        md5hash=md5,
                    )
                    for raw, md5 in PD14_RAWS[url]
                ]
                fragment.validate.validate_spectra(
                    search_path,
                    raw_paths=[i for i in raw_paths],
                    score=20,
                    cpu_count=2,
                    auto_maybe=True,
                )

    def test_load_pd14(self):
        for url, md5 in PD14_URLS:
            with TemporaryDirectory() as tmpdir:
                path = self.fetch_url(
                    PD_URL_BASE + PD14_URL_BASE + url,
                    tmpdir,
                    md5hash=md5,
                )
                search.search.read_search_file(path)

    def test_load_pd22(self):
        for url, md5 in PD22_URLS:
            with TemporaryDirectory() as tmpdir:
                path = self.fetch_url(
                    PD_URL_BASE + PD22_URL_BASE + url,
                    tmpdir,
                    md5hash=md5,
                )
                search.search.read_search_file(path)
