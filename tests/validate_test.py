import hashlib
import logging
import os
import requests
import sys
import tempfile
from unittest import TestCase

try:
    tempfile.TemporaryDirectory
except AttributeError:
    from backports import tempfile

from pycamv import search, main

PD_URL_BASE = (
    "https://media.githubusercontent.com/media/"
    "white-lab/pycamverter-data/master/"
)

PD14_URL_BASE = 'Test%20PD1.4/'
PD22_URL_BASE = 'Test%20PD2.2/'

PD14_URLS = [
    (
        "CK-C1-pY.msf",
        "96eda5b0e47f615cf000d1c5d3ecc8cd",
        [5330, 5732, 5818, 6185, 6296, 7044, 7211, 28124, 28164, 28201, 28257],
    ),
    (
        "FAD-H1-Global.msf",
        "497ba1841faac883619d4f91a86a95cc"
        [1472, 1473, 1518, 1519, 1520, 1521, 1524, 8593, 8594, 8598],
    ),
    (
        "FAD-H1-MPM2.msf",
        "5def08356bcfa5b679835e4a23dd1396",
        [3702, 5739, 5799, 5996, 6088, 6179, 6244, 22278, 23179, 23237, 23332],
    ),
    (
        "Tau-4moHL1-Global.msf",
        "95d8089b5e4657b348bea0868c655478",
        [],
    ),
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
        logging.info('Downloading {} to {}'.format(url, dir))
        response = requests.get(url, stream=True)
        path = os.path.join(dir, url.split('/')[-1])
        path = os.path.splitext(path)[0] + os.path.splitext(path)[1].lower()
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
        for url, md5, include_scans in PD14_URLS:
            with tempfile.TemporaryDirectory() as tmp_dir:
                print(url, md5, tmp_dir)
                logger = logging.getLogger('pycamv')
                logger.setLevel(logging.DEBUG)

                ch = logging.StreamHandler(sys.stdout)
                ch.setLevel(logging.DEBUG)

                search_path = self.fetch_url(
                    PD_URL_BASE + PD14_URL_BASE + url,
                    tmp_dir,
                    md5hash=md5,
                )
                raw_paths = [
                    self.fetch_url(
                        PD_URL_BASE + PD14_URL_BASE + raw,
                        tmp_dir,
                        md5hash=md5,
                    )
                    for raw, md5 in PD14_RAWS[url]
                ]
                main.main(
                    [
                        '--score=20',
                        '--cpus=1',
                        '--scans',
                        include_scans,
                        search_path,
                    ] + [i for i in raw_paths],
                )

    def test_load_pd14(self):
        for url, md5 in PD14_URLS:
            with tempfile.TemporaryDirectory() as tmp_dir:
                path = self.fetch_url(
                    PD_URL_BASE + PD14_URL_BASE + url,
                    tmp_dir,
                    md5hash=md5,
                )
                search.search.read_search_file(path)

    def test_load_pd22(self):
        for url, md5 in PD22_URLS:
            with tempfile.TemporaryDirectory() as tmp_dir:
                path = self.fetch_url(
                    PD_URL_BASE + PD22_URL_BASE + url,
                    tmp_dir,
                    md5hash=md5,
                )
                search.search.read_search_file(path)
