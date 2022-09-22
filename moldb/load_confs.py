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
Loads 3D conformers from a file into the enumeration table.
The input data should be generated with the conformers.py script.
It will be in the form of a tab delineated file that has a multiple 3D conformers for each enumerated state.
The file should contain these fields:

* the 3D coordinates section from the cxmiles
* enumerated ID
* energy
* energy delta (difference to lowest energy conformer)

The coordinates bit must be the coordiantes section of the cxsmiles extension. e.g if the cxsmiles looked like this:
  CCO |(2.84126,-0.0586636,2.68552;3.04156,1.20289,1.86426;3.754,0.913667,0.663692)|
then the first field must be like this:
  2.84126,-0.0586636,2.68552;3.04156,1.20289,1.86426;3.754,0.913667,0.663692

Generate data like this using the conformers.py module with the -coords-only option.
"""

import argparse

from sqlalchemy import delete
from sqlalchemy.orm import Session

from . import models
from . models import Enumeration, Conformer

from dm_job_utilities.dm_log import DmLog


engine = models.get_engine(echo=False)


def _load_confs(session, moldata):
    items = []
    for data in moldata:
        id = data[0]
        coords = data[1]
        energy = float(data[2])
        energy_delta = float(data[3])
        print('loading', data)
        items.append(Conformer(enumeration_id=id, coords=coords, energy=energy, energy_delta=energy_delta))

    session.bulk_save_objects(items)


def load_data(inputs, chunk_size=100, interval=None):
    """
    Loads enumerated 3D molecules from a SD-file into the enumeration table.

    :param input: file name to load (contents: cords\tid)
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
                    coords = t[0].strip()
                    id = t[1].strip()
                    energy = t[3].strip()
                    energy_delta = t[4].strip()

                    count += 1
                    recordno += 1

                    # if it's a new id then wipe all existing rows with that id from the db
                    if id != current_id:
                        current_id = id
                        num_ids += 1
                        stmt = delete(Conformer).where(Conformer.enumeration_id == id)
                        session.execute(stmt)

                    if interval and recordno % interval == 0:
                        DmLog.emit_event("Processed {} records".format(recordno))

                    chunk.append((id, coords, energy, energy_delta))
                    # if chunk has got to the right size then do a bulk insert into the db
                    if len(chunk) == chunk_size:
                        _load_confs(session, chunk)
                        chunk = []

                # load any dangling chunk
                if len(chunk) > 0:
                    _load_confs(session, chunk)
                    chunk = []

                    session.commit()

    DmLog.emit_event("Total of {} records loaded, {} distinct IDs".format(count, num_ids))
    return count


def main():

    # Example:
    #   python -m moldb.load_enums -i enumerations.sdf --interval 10000
    #   python -m moldb.load_enums -i enumerations.cxsmi --interval 10000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='load_enums')
    parser.add_argument('-i', '--input', nargs='+', help="Input file (.cxsmi")
    parser.add_argument("--chunk-size", type=int, default=100, help="Bulk insert chunk size")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("load_confs: ", args)

    load_data(args.input, interval=args.interval)


if __name__ == "__main__":
    main()
