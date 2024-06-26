# Copyright 2022 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Loads enumerated 3D molecules from a file of cxsmiles into the enumeration table.
The input data should be generated with the enumerate.py script.
It will be in the form of a tab delineated file that has a single 3D conformer of each enumerated state.
The file should contain these fields:

* cxsmiles (with 3D coordinates)
* molecule ID
* enumeration code (B=base structure, T=tautomer, M=microstate, C=stereoisomer)

For each molecule ID all enumerated molecules are first deleted from the DB before
the new data is loaded. This means it is important to be consistent about which forms are being
enumerated or the db will contain inconsistent enumerated forms. Unless you have particular
reasons otherwise, you should use enumerate.py with the --enumerate-charges, --enumerate-chirals
and --enumerate-tautomers options but not the --combinatorial option.

For a particular db apply the same enumeration settings to all uses of enumerate.py.
"""

import argparse

from sqlalchemy import delete
from sqlalchemy.orm import Session

from . import models
from . models import Enumeration

from dm_job_utilities.dm_log import DmLog


engine = models.get_engine(echo=False)


def _load_mols(session, moldata):
    items = []
    for data in moldata:
        items.append(Enumeration(molecule_id=data[0], smiles=data[1], coords=data[2], code=data[3]))

    session.bulk_save_objects(items)


def load_data(inputs, chunk_size=100, purge=False, interval=None):
    """
    Loads enumerated 3D molecules from a SD-file into the enumeration table.

    :param input: file name to load (contents: cxsmi\tid\tcode
    :param chunk_size: Bulk insert chunk size
    :param interval: Reporting interval
    :return:
    """

    count = 0
    chunk = []
    current_id = None
    num_ids = 0

    with Session(engine) as session:
        for input in inputs:
            recordno = 0
            DmLog.emit_event("Handling", input)

            with open(input) as f:
                for line in f:
                    if not line:
                        break
                    t = line.split('\t')
                    cxsmi = t[0].strip()
                    id = t[1].strip()
                    code = t[2].strip()

                    parts = cxsmi.split(' ')
                    smi = parts[0]
                    coords = parts[1][1:-1]

                    count += 1
                    recordno += 1

                    if id != current_id:
                        current_id = id
                        num_ids += 1
                        # if purge then wipe all existing rows with that id from the db
                        if purge:
                            stmt = delete(Enumeration).where(Enumeration.molecule_id == id)
                            session.execute(stmt)

                    if interval and recordno % interval == 0:
                        DmLog.emit_event("Processed {} records".format(recordno))

                    chunk.append((id, smi, coords, code))
                    # if chunk has got to the right size then do a bulk insert into the db
                    if len(chunk) == chunk_size:
                        _load_mols(session, chunk)
                        chunk = []

                # load any dangling chunk
                if len(chunk) > 0:
                    _load_mols(session, chunk)
                    chunk = []

                    session.commit()

    DmLog.emit_event("Total of {} records loaded, {} distinct IDs".format(count, num_ids))
    return count


def main():

    # Example:
    #   python -m moldb.load_enums -i enumerations.cxsmi --interval 10000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='load_enums')
    parser.add_argument('-i', '--input', nargs='+', help='Input file (.cxsmi)')
    parser.add_argument('-c', '--chunk-size', type=int, default=100, help='Bulk insert chunk size')
    parser.add_argument('-p', '--purge', action='store_true', help='Purge data before loading')
    parser.add_argument('--interval', type=int, help='Reporting interval')

    args = parser.parse_args()
    DmLog.emit_event("load_enums: ", args)

    load_data(args.input, purge=args.purge, interval=args.interval)


if __name__ == "__main__":
    main()
