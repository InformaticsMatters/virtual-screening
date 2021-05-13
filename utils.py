import sys
from hashlib import sha256


default_num_chars = 2
default_num_levels = 2


def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


def get_path_from_digest(digest, num_chars=default_num_chars, num_levels=default_num_levels):
    parts = []
    start = 0
    for l in range(0, num_levels):
        end = start + num_chars
        p = digest[start:end]
        parts.append(p)
        start = start + num_chars
    return parts
    
