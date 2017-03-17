# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys

from pycamv import validate, export, gui, scan_list, search, __version__


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
        '-V', '--version',
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "--show_gui",
        action="store_true",
        help="Show GUI for converting files.",
    )
    parser.add_argument(
        "--raw_path",
    )
    parser.add_argument(
        "--search_path",
        help="MASCOT or ProteomeDiscoverer search files."
    )
    parser.add_argument(
        "--scans_path",
    )
    parser.add_argument(
        "--scans",
        nargs="*", type=int,
        help="Individual scans to select for validation.",
    )
    parser.add_argument(
        "--out_path",
    )
    parser.add_argument(
        'files', nargs='*',
    )
    return parser.parse_args(args)


def main(args):
    args = _parse_args(args)

    if not args.verbose:
        level = logging.WARNING
    elif args.verbose in [1]:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG

    logging.basicConfig(
        # filename="pycamv.log",
        level=level,
    )

    LOGGER.debug(sys.argv)
    LOGGER.debug(args)

    if args.show_gui:
        gui.run_gui(args)
    else:
        if args.search_path is None:
            searches = [
                i for i in args.files
                if os.path.splitext(i.lower())[1] in SEARCH_EXTS
            ]

            if len(searches) == 1:
                args.search_path = searches[0]

        if args.raw_path is None:
            raws = [
                i for i in args.files
                if os.path.splitext(i.lower())[1] in RAW_EXTS
            ]

            if len(raws) == 1:
                args.raw_path = raws[0]

        if args.scans_path is None:
            scans = [
                i for i in args.files
                if os.path.splitext(i.lower())[1] in SCANS_EXTS
            ]

            if len(scans) == 1:
                args.scans_path = scans[0]

        if (
            args.search_path is None or
            args.raw_path is None
        ):
            raise Exception(
                "Missing either input search / raw paths"
            )

        options, peak_hits, scan_mapping, precursor_windows, label_windows = (
            validate.validate_spectra(
                search_path=args.search_path,
                raw_path=args.raw_path,
                scans_path=args.scans_path,
                scan_list=args.scans,
            )
        )

        if args.out_path is None:
            args.out_path = os.path.splitext(args.raw_path)[0] + ".camv.gz"

        export.export_to_camv(
            args.out_path,
            peak_hits, scan_mapping, precursor_windows, label_windows,
        )

if __name__ == "__main__":
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
