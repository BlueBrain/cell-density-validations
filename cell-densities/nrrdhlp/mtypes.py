import re

mtype_regex = r"L\d+_[A-Z:]{2,6}"


def parse_mtype(fn):
    for match in re.findall(mtype_regex, fn):
        return match
    return None
