
import re

RE_DISCOVERER_DESCRIPTION = re.compile(
    r"^>sp\|([\dA-Za-z]+)\|([\dA-Za-z_]+) (.+) OS=(.+)$"
)
RE_MASCOT_DESCRIPTION = re.compile(
    r"^(.+) OS=(.+)$"
)
RE_DYN_MODS = re.compile(r"((\d+) )?(.+) \((.+)\)")
RE_BY_ION_POS = re.compile(r"([abcxyz])_\{(\d+)\}")
RE_B_Y_IONS = re.compile(r"([abcxyz]_\{[0-9]+\})(.*)\^\{\+\}")
RE_SCAN_NUM = re.compile(r"(scans:|Cmpd_)(\d+)")
RE_COLLISION_TYPE = re.compile(r".*@([A-Za-z]+)\d+")
RE_PRECURSOR_SCAN = re.compile(r"scan=(\d+)")
