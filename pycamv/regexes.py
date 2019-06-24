
import re

RE_DISCOVERER_DESCRIPTION = re.compile(
    r"^>([A-Za-z]+)\|([\dA-Za-z]+)\|([\dA-Za-z_]+) (.+) OS=(.+)$"
)
"""
Regex matching Discoverer protein description.
"""
RE_FALLBACK_DESCRIPTION = re.compile(
    r"^>([A-Za-z]+)\|([\dA-Za-z]+)\|([\dA-Za-z_]+)\|([\dA-Za-z_\.]+)\|([^|]+?)( \[(.+)\])?$"
)
"""
Regex matching alternative protein description.
"""
RE_MASCOT_DESCRIPTION = re.compile(
    r"^(.+) OS=(.+)$"
)
"""
Regex matching MASCOT protein description.
"""
RE_DYN_MODS = re.compile(r"((\d+) )?(.+) \((.+)\)")
"""
Regex matching dynamic modifications.
"""
RE_BY_ION_POS = re.compile(r"([abcxyz])_\{(\d+)\}")
"""
Regex matching b/y ion position.
"""
RE_B_Y_IONS = re.compile(r"([abcxyz]_\{[0-9]+\})(.*)\^\{\+\}")
"""
Regex matching b/y ion names.
"""
RE_SCAN_NUM = re.compile(r"(scans:|Cmpd_)(\d+)")
"""
Regex matching Scan Number in .mzML.
"""
RE_COLLISION_TYPE = re.compile(r".*@([A-Za-z]+)\d+")
"""
Regex matching Collision Type in .mzML.
"""
RE_PRECURSOR_SCAN = re.compile(r"scan=(\d+)")
"""
Regex matching Precursor Scan name in .mzML.
"""
