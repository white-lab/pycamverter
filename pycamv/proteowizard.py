"""
Provides functions for interacting with MS data through ProteoWizard.
"""

from __future__ import absolute_import, division

import cgi
import logging
import os
import platform
import shutil
import subprocess
import sys
import tempfile

import requests
import pymzml


LOGGER = logging.getLogger("pycamv.proteowizard")


THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROTEOWIZARD_DIR = os.path.join(THIS_DIR, "ProteoWizard")

PROTEOWIZARD_VERSION = "3.0.10577"
PROTEOWIZARD_PATH = os.path.join(
    PROTEOWIZARD_DIR,
    "ProteoWizard {}".format(PROTEOWIZARD_VERSION),
)

PROTEOWIZARD_MSI_32_URLS = [
    (
        (
            "https://www.dropbox.com/s/we46o1hvkzlpel7/"
            "pwiz-setup-3.0.10577-x86.msi?dl=1"
        ),
        "GET",
        None,
    ),
    (
        "http://data.mallicklab.com/download.php",
        "POST",
        {"downloadtype": "bt36i"},
    ),
]

PROTEOWIZARD_MSI_64_URLS = [
    (
        (
            "https://www.dropbox.com/s/9rawc44u84jr9ip/"
            "pwiz-setup-3.0.10577-x86_64.msi?dl=1"
        ),
        "GET",
        None,
    ),
    (
        "http://data.mallicklab.com/download.php",
        "POST",
        {"downloadtype": "bt83i"},
    ),
]


def fetch_proteowizard(urls=None):
    LOGGER.debug("Proteowizard: {}".format(PROTEOWIZARD_PATH))

    if os.path.exists(PROTEOWIZARD_PATH):
        return

    LOGGER.info("ProteoWizard not installed, fetching now.")

    if platform.system() not in ["Windows"]:
        raise Exception("Proteowizard install not supported on your platform")

    if urls is None:
        if platform.architecture()[0] == "64bit":
            urls = PROTEOWIZARD_MSI_64_URLS
        else:
            urls = PROTEOWIZARD_MSI_32_URLS

    tmpdir = tempfile.mkdtemp()

    # Download the .msi file
    responses = []
    for url, req_type, post in urls:
        if req_type in ["GET"]:
            response = requests.get(url, stream=True)
        elif req_type in ["POST"]:
            response = requests.post(url, data=post, stream=True)
        else:
            raise Exception("Unknown request type: {}".format(req_type))

        if not response.ok:
            responses.append(response)
            continue

        out_path = os.path.join(
            tmpdir,
            cgi.parse_header(
                response.headers["Content-Disposition"]
            )[-1]["filename"]
        )

        with open(out_path, mode="wb") as f:
            for block in response.iter_content(1024):
                # print(block)
                # return
                f.write(block)

            break
    else:
        raise Exception("Unable to download file: {}".format(responses))

    # Extract the msi file's contents
    extract_path = os.path.join(tmpdir, "msi_extract")
    cmd = [
        "msiexec",
        "/a",
        out_path,
        "/qb",
        "TARGETDIR=\"{}\"".format(extract_path),
    ]
    subprocess.check_call(" ".join(cmd), shell=True)

    # Copy the msi file's contents to PROTEOWIZARD_DIR
    src = os.path.join(
        extract_path,
        "PFiles",
        "ProteoWizard",
    )

    shutil.rmtree(PROTEOWIZARD_DIR, ignore_errors=True)
    shutil.copytree(src, PROTEOWIZARD_DIR)
    shutil.rmtree(tmpdir)


def raw_to_mzml(raw_path, out_dir, scans=None, mz_window=None):
    """
    Covert a RAW file to .mzML using ProteoWizard.

    Parameters
    ----------
    raw_path : str
    out_dir : str
    scans : list of int, optional
    mz_window : list of int, optional

    Returns
    -------
    :class:`pymzml.run.Reader<run.Reader>`
    """
    fetch_proteowizard()

    ms_convert_path = os.path.join(PROTEOWIZARD_PATH, "msconvert.exe")

    # Create a config file,
    config, config_path = tempfile.mkstemp(suffix=".txt", text=True)

    with os.fdopen(config, "w+") as config:
        if scans:
            config.write(
                "filter=\"scanNumber {}\"\n".format(
                    " ".join(str(scan) for scan in scans)
                )
            )

        if mz_window:
            config.write(
                "filter=\"mzWindow [{},{}]\"\n".format(
                    mz_window[0], mz_window[1],
                )
            )

    # Run msconvert to convert raw file to mzML
    LOGGER.info("Converting \"{}\" to .mzML format.".format(raw_path))

    cmd = [
        ms_convert_path,
        raw_path,
        "-o", out_dir,
        "--mzML",
        "-c", config_path,
    ]

    out = subprocess.check_output(cmd)
    LOGGER.debug(out.decode(sys.stdout.encoding))

    os.remove(config_path)

    # Read the file into memory using pymzml
    basename = os.path.splitext(os.path.basename(raw_path))[0]
    out_path = os.path.join(out_dir, "{}.mzML".format(basename))
    data = pymzml.run.Reader(
        out_path,
        extraAccessions=[
            ("MS:1000827", ["value"]),  # isolation window target m/z
            ("MS:1000828", ["value"]),  # isolation window lower offset
            ("MS:1000829", ["value"]),  # isolation window upper offset
            ("MS:1000512", ["value"]),  # filter string
        ],
    )

    return data
