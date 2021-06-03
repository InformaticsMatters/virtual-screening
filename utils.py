import io
from datetime import datetime, timezone
import logging
import sys

default_num_chars = 2
default_num_levels = 2

_SBUF = io.StringIO()
_INFO = logging.getLevelName(logging.INFO)

def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


def log_dm_event(*args):
    """Generate a Data Manager-compliant event message.
    The Data Manager watches stdout and interprets certain formats
    as an 'event'. These are then made available to the client.
    Here we write the message using the expected format.
    """
    _ = _SBUF.truncate(0)
    print(*args, file=_SBUF)
    msg_time = datetime.now(timezone.utc).replace(microsecond=0)
    print('%s # %s -EVENT- %s' % (msg_time.isoformat(),
                                  _INFO,
                                  _SBUF.getvalue().strip()))


def get_path_from_digest(digest, num_chars=default_num_chars, num_levels=default_num_levels):
    parts = []
    start = 0
    for l in range(0, num_levels):
        end = start + num_chars
        p = digest[start:end]
        parts.append(p)
        start = start + num_chars
    return parts
    
