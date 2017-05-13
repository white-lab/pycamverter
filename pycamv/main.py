# -*- coding: utf-8 -*-

import argparse
import logging
import multiprocessing
import os
import sys

from pycamv import validate, scan_list, search, __version__


LOGGER = logging.getLogger("pycamv.main")

RAW_EXTS = [".raw", ".mgf", ".d", ".wiff"]
SEARCH_EXTS = list(search.BACKENDS.keys())
SCANS_EXTS = list(scan_list.BACKENDS.keys())


def _parse_args(args):
    """
    Parses arguments from a argv format.

    Parameters
    ----------
    args : list of str

    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="PyCAMVerter",
        description="Aww yeah, mass specs!",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="count",
        help="Increase verbosity of output.",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="count",
        help="Decrease verbosity of output.",
    )
    parser.add_argument(
        '-V', '--version',
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "--cpus",
        type=int,
        help="Limit the number of concurrent processes."
    )
    parser.add_argument(
        "--reprocess",
        action="store_true",
        help="Reprocess a set of scans, without limiting ptm combinations.",
    )
    parser.add_argument(
        "--no-auto-maybe",
        action="store_false",
        dest="auto_maybe",
        help="Don't auto-assign all peptides with best MASCOT rank as 'maybe'",
    )
    parser.add_argument(
        "--raw-paths",
        nargs="+",
        help="Raw data file(s) containing mass spec data.",
    )
    parser.add_argument(
        "--search-path",
        help="MASCOT or ProteomeDiscoverer search files."
    )
    parser.add_argument(
        "--scans-path",
        help=".xlsx or .csv file listing scans to select for validation.",
    )
    parser.add_argument(
        "--scans",
        nargs="*", type=int,
        help="Individual scans to select for validation.",
    )
    parser.add_argument(
        "--mat-sessions",
        nargs="+",
        help="Path to CAMV-Matlab session files.",
    )
    parser.add_argument(
        "--out-path",
        help="Output path for CAMV export.",
    )
    parser.add_argument(
        'files',
        nargs='*',
        help="Raw, search, or scan list files, determined by file extension."
    )
    return parser, parser.parse_args(args)


def main(args):
    parser, args = _parse_args(args)

    verbosity = (args.verbose or 0) - (args.quiet or 0)

    if verbosity <= -2:
        level = logging.CRITICAL
    elif verbosity == -2:
        level = logging.ERROR
    elif verbosity == -1:
        level = logging.WARNING
    elif verbosity == 0:
        level = logging.INFO
    elif verbosity > 0:
        level = logging.DEBUG

    logger = logging.getLogger('pycamv')
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    fh = logging.FileHandler('pycamv.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    LOGGER.debug(sys.argv)
    LOGGER.debug(args)

    if args.search_path is None:
        searches = [
            i for i in args.files
            if os.path.splitext(i.lower())[1] in SEARCH_EXTS
        ]

        if len(searches) == 1:
            args.search_path = searches[0]

    if args.raw_paths is None:
        raws = [
            i for i in args.files
            if os.path.splitext(i.lower())[1] in RAW_EXTS
        ]

        if len(raws) > 0:
            if args.raw_paths:
                args.raw_paths += raws
            else:
                args.raw_paths = raws

    if args.scans_path is None:
        scans = [
            i for i in args.files
            if os.path.splitext(i.lower())[1] in SCANS_EXTS
        ]

        if len(scans) == 1:
            args.scans_path = scans[0]

    if (
        args.search_path is None
    ):
        parser.print_help()
        raise Exception(
            "Missing search input path"
        )

    validate.validate_spectra(
        search_path=args.search_path,
        raw_paths=args.raw_paths,
        scans_path=args.scans_path,
        scan_list=args.scans,
        mat_sessions=args.mat_sessions,
        out_path=args.out_path,
        reprocess=args.reprocess,
        auto_maybe=args.auto_maybe,
        cpu_count=args.cpus,
    )


if __name__ == "__main__":
    multiprocessing.freeze_support()
    # is_frozen_executable = getattr(sys, u'frozen', False)
    #
    # if is_frozen_executable:
    #     os.chdir(os.path.dirname(sys.executable))
    # else:
    #     os.chdir(os.path.join(os.path.dirname(__file__), u'..'))

    try:
        main(sys.argv[1:])
    except Exception as e:
        LOGGER.error("PyCAMV Converter has crashed!", exc_info=True)
        raise e
