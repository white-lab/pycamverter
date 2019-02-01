"""
Provides functions for interacting with MS data through ProteoWizard.
"""

from __future__ import absolute_import, division

import hashlib
import logging
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import requests

try:
    tempfile.TemporaryDirectory
except AttributeError:
    from backports import tempfile

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
                    "https://www.dropbox.com/s/j54gb26dpsja6tn/"
                    "pwiz-bin-windows-x86_64-vc141-release-3_0_19014_f9d5b8a3b"
                    ".tar.bz2?dl=1"
                ),
                "GET",
                None,
            ),
            (
                "http://teamcity.labkey.org/guestAuth/app/rest/builds/"
                "id:686855/artifacts/content/pwiz-bin-windows-x86_64"
                "-vc141-release-3_0_19014_f9d5b8a3b.tar.bz2",
                "GET",
                None,
            ),
        ],
        "5b1af1f0817a45f9600bb459993abb2b",
    ),
    "win64": (
        [
            (
                (
                    "https://www.dropbox.com/s/hsd3ew8lm8gbk23/"
                    "pwiz-bin-windows-x86-vc141-release-3_0_19014_f9d5b8a3b"
                    ".tar.bz2?dl=1"
                ),
                "GET",
                None,
            ),
            (
                "http://teamcity.labkey.org/guestAuth/app/rest/builds/"
                "id:686846/artifacts/content/pwiz-bin-windows-x86"
                "-vc141-release-3_0_19014_f9d5b8a3b.tar.bz2",
                "GET",
                None,
            ),
        ],
        "247e28638d498050965482c8c0faa4ea"
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
    LOGGER.debug("Proteowizard: {}".format(PROTEOWIZARD_DIR))

    if os.path.exists(PROTEOWIZARD_DIR):
        return

    LOGGER.info(
        "ProteoWizard not installed, fetching to {}.".format(PROTEOWIZARD_DIR)
    )

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

    shutil.rmtree(PROTEOWIZARD_DIR, ignore_errors=True)

    # Copy the msi file's contents to PROTEOWIZARD_DIR
    with tarfile.open(out_path, "r:bz2") as f:
        f.extractall(PROTEOWIZARD_DIR)

    shutil.rmtree(tmpdir)


def _write_config(config_path, scans=None, mz_window=None):
    with open(config_path, "w+") as config:
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


def raw_to_mzml(raw_path, scans=None, mz_window=None):
    """
    Covert a RAW file to .mzML using ProteoWizard.

    Parameters
    ----------
    raw_path : str
    scans : list of int, optional
    mz_window : list of int, optional

    Returns
    -------
    :class:`pymzml.run.Reader<run.Reader>`
    """
    LOGGER.info("Converting \"{}\" to .mzML format.".format(raw_path))

    basename = os.path.splitext(os.path.basename(raw_path))[0]
    config_name = '{}_msconvert.txt'.format(
        os.path.splitext(os.path.basename(raw_path))[0]
    )
    tmp_dir = None

    if platform.system() in ["Windows"]:
        fetch_proteowizard()

        cmd = [
            os.path.join(PROTEOWIZARD_DIR, "msconvert.exe")
        ]

        tmp_dir = tempfile.TemporaryDirectory()

        out_dir = tmp_dir.name
        config_dir = tmp_dir.name
        mzml_path = os.path.join(tmp_dir.name, "{}.mzML".format(basename))
    else:
        raw_dir = os.path.dirname(raw_path)
        cmd = [
            'docker',
            'run',
            '-t',
            '-v', '{}:/data'.format(raw_dir),
            'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:x64',
            'wine',
            'msconvert',
        ]

        out_dir = '/data'
        config_dir = raw_dir
        mzml_path = os.path.join(raw_dir, "{}.mzML".format(basename))
        raw_path = os.path.join(out_dir, os.path.basename(raw_path))

    _write_config(
        os.path.join(config_dir, config_name),
        scans=scans,
        mz_window=mz_window,
    )

    # Run msconvert to convert raw file to mzML
    cmd += [
        raw_path,
        "-o", out_dir,
        "--mzML",
        "-c", os.path.join(out_dir, config_name),
    ]

    encoding = sys.stdout.encoding or "utf-8"

    LOGGER.debug('Calling subprocess: {}'.format(" ".join(cmd)))

    try:
        out = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as err:
        LOGGER.error(
            "Error Running msconvert:\n{}".format(err.output.decode(encoding))
        )
        raise

    LOGGER.debug(out.decode(encoding))

    # Read the file into memory using pymzml
    data = pymzml.run.Reader(
        mzml_path,
        extraAccessions=[
            ("MS:1000827", ["value"]),  # isolation window target m/z
            ("MS:1000828", ["value"]),  # isolation window lower offset
            ("MS:1000829", ["value"]),  # isolation window upper offset
            ("MS:1000512", ["value"]),  # filter string
        ],
    )
    data._tmp_dir = tmp_dir

    return data
