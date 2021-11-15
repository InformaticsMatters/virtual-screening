import io
from datetime import datetime, timezone
import logging
import numbers
import sys, os

default_num_chars = 2
default_num_levels = 2

_SBUF = io.StringIO()
_INFO = logging.getLevelName(logging.INFO)

# A value used to ensure cost lines are unique (regardless of timestamp)
# It's incremented on entry to log_dm_cost()
_DM_COST_SEQUENCE_NUMBER = 0

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

def log_dm_cost(cost, cumulative=True):
    """Generate a Data Manager-compliant cost message.
    The Data Manager watches stdout and interprets certain formats
    as a 'cost', typically used for billing purposes.

    The cost must be a non-negative number. It is assumed to be cumulative
    unless cumulative is set to False. Cumulative values are written
    with a '+' prefix, i.e. '+1' or '+0'. Non-cumulative (absolute) costs
    are written without a '+' prefix, i.e. '1' or '0'.
    """
    global _DM_COST_SEQUENCE_NUMBER

    # Ensure this cost message is unique
    _DM_COST_SEQUENCE_NUMBER += 1

    # Cost is always expected to be a number that's not negative.
    assert isinstance(cost, numbers.Number)
    assert cost >= 0

    cost_str = str(cost)
    if cumulative:
        cost_str = '+' + cost_str
    msg_time = datetime.now(timezone.utc).replace(microsecond=0)
    print('%s # %s -COST- %s %d' % (msg_time.isoformat(),
                                    _INFO,
                                    cost_str,
                                    _DM_COST_SEQUENCE_NUMBER))

def get_path_from_digest(digest, num_chars=default_num_chars, num_levels=default_num_levels):
    parts = []
    start = 0
    for l in range(0, num_levels):
        end = start + num_chars
        p = digest[start:end]
        parts.append(p)
        start = start + num_chars
    return parts


def expand_path(path):
    """
    Create any necessary directories to ensure that the file path is valid
    
    :param path: a filename or directory that might or not exist
    """
    head_tail = os.path.split(path)
    if head_tail[0]:
        if not os.path.isdir(head_tail[0]):
            log('Creating directories for', head_tail[0])
            os.makedirs(head_tail[0], exist_ok=True)

def UpdateChargeFlagInAtomBlock(mb):
    """
    See https://sourceforge.net/p/rdkit/mailman/message/36425493/
    """
    f="{:>10s}"*3+"{:>2}{:>4s}"+"{:>3s}"*11
    chgs = []    # list of charges
    lines = mb.split("\n")
    if mb[0] == '' or mb[0] == "\n":
        del lines[0]
    CTAB = lines[2]
    atomCount = int(CTAB.split()[0])
    # parse mb line per line
    for l in lines:
        # look for M CHG property
        if l[0:6] == "M  CHG":
            records = l.split()[3:]    # M  CHG X is not needed for parsing, the info we want comes afterwards
            # record each charge into a list
            for i in range(0,len(records),2):
                idx = records[i]
                chg = records[i+1]
                chgs.append((int(idx), int(chg)))    # sort tuples by first element?
            break    # stop iterating

    # sort by idx in order to parse the molblock only once more
    chgs = sorted(chgs, key=lambda x: x[0])

    # that we have a list for the current molblock, attribute each charges
    for chg in chgs:
        i=3
        while i < 3+atomCount:    # do not read from beginning each time, rather continue parsing mb!
            # when finding the idx of the atom we want to update, extract all fields and rewrite whole sequence
            if i-2 == chg[0]:    # -4 to take into account the CTAB headers, +1 because idx begin at 1 and not 0
                fields = lines[i].split()
                x=fields[0]
                y=fields[1]
                z=fields[2]
                symb=fields[3]
                massDiff=fields[4]
                charge=fields[5]
                sp=fields[6]
                hc=fields[7]
                scb=fields[8]
                v=fields[9]
                hd=fields[10]
                nu1=fields[11]
                nu2=fields[12]
                aamn=fields[13]
                irf=fields[14]
                ecf=fields[15]
                # update charge flag
                if chg[1] == -1:
                    charge = '5'
                elif chg[1] == -2:
                    charge = '6'
                elif chg[1] == -3:
                    charge = '7'
                elif chg[1] == 1:
                    charge = '3'
                elif chg[1] == 2:
                    charge = '2'
                elif chg[1] == 3:
                    charge = '1'
                else:
                    print("ERROR! " + str(lines[0]) + "unknown charge flag: " + str(chg[1]))    # print name then go to next chg
                    break
                # update modatom block line
                lines[i] = f.format(x,y,z,symb,massDiff,charge,sp,hc,scb,v,hd,nu1,nu2,aamn,irf,ecf)
            i+=1
    #print("\n".join(lines))
    del lines[-1]    # remove empty element left because last character before $$$$ is \n
    upmb = "\n" + "\n".join(lines)
    return(upmb)


def read_delimiter(input):
    if input:
        if 'tab' == input:
            delimiter = '\t'
        elif 'space' == input:
            delimiter = None
        elif 'comma' == input:
            delimiter = ','
        elif 'pipe' == input:
            delimiter = '|'
        else:
            delimiter = input
    else:
        delimiter = None
    return delimiter
