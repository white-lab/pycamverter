"""
Provides functions for interacting with MS data through ProteoWizard.
"""

from __future__ import absolute_import, division

import cgi
import hashlib
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

PROTEOWIZARD_VERSION = "3.0.19014"
PROTEOWIZARD_PATH = os.path.join(
    PROTEOWIZARD_DIR,
    "ProteoWizard {}".format(PROTEOWIZARD_VERSION),
)

# URLs + MD5 Hash:
# Windows Powershell: CertUtil -hashfile <filename> MD5
# Linux Commandline: md5sum <filename>
PROTEOWIZARD_MSI_URLS = {
    "win32": (
        [
            (
                (
                    "https://www.dropbox.com/s/43jfkvqu268gc55/"
                    "pwiz-setup-3.0.19014.f9d5b8a3b-x86_64.msi?dl=1"
                ),
                "GET",
                None,
            ),
            (
                "http://teamcity.labkey.org/guestAuth/app/rest/builds/"
                "id:686855/artifacts/content/"
                "pwiz-setup-3.0.19014.f9d5b8a3b-x86_64.msi",
                "GET",
                None,
            ),
        ],
        "5b812f1cc1395e37819a53412e8e8c79",
    ),
    "win64": (
        [
            (
                (
                    "https://www.dropbox.com/s/39mkw5wjcjf8ft3/"
                    "pwiz-setup-3.0.19014.f9d5b8a3b-x86.msi?dl=1"
                ),
                "GET",
                None,
            ),
            (
                "http://teamcity.labkey.org/guestAuth/app/rest/builds/"
                "id:686846/artifacts/content/"
                "pwiz-setup-3.0.19014.f9d5b8a3b-x86.msi",
                "GET",
                None,
            ),
        ],
        "d2f1ecca398bb7f10ac21815ac63b3d9"
    )
}


def fetch_proteowizard(urls=None, md5hash=None):
    """
    Download ProteoWizard to this module's directory.

    Parameters
    ----------
    urls : list of str, optional
        URL for ProteoWizard installer.
    md5hash : list of str, optional
        MD5 Hash for ProteoWizard installer.
    """
    LOGGER.debug("Proteowizard: {}".format(PROTEOWIZARD_PATH))

    if os.path.exists(PROTEOWIZARD_PATH):
        return

    LOGGER.info("ProteoWizard not installed, fetching now.")

    if platform.system() not in ["Windows"]:
        raise Exception("Proteowizard install not supported on your platform")

    if urls is None:
        urls, md5hash = PROTEOWIZARD_MSI_URLS[
            "win64"
            if platform.architecture()[0] == "64bit" else
            "win32"
        ]

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

        try:
            filename = cgi.parse_header(
                response.headers["Content-Disposition"]
            )[-1]["filename"]
        except KeyError:
            filename = url.rsplit("/", 1)[1].rsplit("?")[0]

        out_path = os.path.join(
            tmpdir,
            filename,
        )
        hash_md5 = hashlib.md5()

        with open(out_path, mode="wb") as f:
            for block in response.iter_content(1024):
                hash_md5.update(block)
                f.write(block)

        if md5hash is not None and hash_md5.hexdigest() != md5hash:
            LOGGER.warning(
                "MD5 hash of {} does not match record: {} != {}"
                .format(url, md5hash, hash_md5.hexdigest())
            )
            os.remove(out_path)
            continue

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

    try:
        out = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as err:
        LOGGER.error("Error Running msconvert:\n{}".format(err.output))
        raise

    encoding = sys.stdout.encoding or "utf-8"
    LOGGER.debug(out.decode(encoding))

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
