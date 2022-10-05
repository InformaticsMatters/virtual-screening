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

import argparse, time, os, random, string, traceback

import utils
from . import models

from sqlalchemy import text
from sqlalchemy.orm import Session
from sqlalchemy import MetaData, Table, Column, Integer, SmallInteger, String, Text, Float

from dm_job_utilities.dm_log import DmLog

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen

engine = models.get_engine(echo=False)

def extract_mols(outfile, count=None):

    utils.expand_path(outfile)
    DmLog.emit_event('Writing to', outfile)

    with Session(engine) as session:

        if count:
            sql = "COPY (SELECT smiles, id FROM molecule WHERE hac IS NULL LIMIT {}) TO STDOUT".format(count)
        else:
            sql = "COPY (SELECT smiles, id FROM molecule WHERE hac IS NULL) TO STDOUT"

        DmLog.emit_event('Extracting uncalculated molecules using:\n', sql)

        output = open(outfile, 'w')
        cursor = session.connection().connection.cursor()
        cursor.copy_expert(sql, output)

        session.commit()


def main():

    ### typical usage:
    # python -m moldb.extract_need_molprops --outfile foo.smi --count 1000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Extract molecules needing property calculation')
    parser.add_argument("--outfile", required=True, help="Output file name")
    parser.add_argument("--count", type=int, help="Max number of molecules to work on")

    args = parser.parse_args()
    DmLog.emit_event("extract_need_molprops: ", args)

    t0 = time.time()
    extract_mols(args.outfile, count=args.count)
    t1 = time.time()
    # Duration? No less than 1 second?
    duration_s = int(t1 - t0)
    if duration_s < 1:
        duration_s = 1

    DmLog.emit_event('Processing complete in {} seconds'.format(duration_s))
    # DmLog.emit_cost(count)


if __name__ == "__main__":
    main()
