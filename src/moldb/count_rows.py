#!/usr/bin/env python

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

import argparse

from sqlalchemy import text

from . import models

from dm_job_utilities.dm_log import DmLog


engine = models.get_engine(echo=False)


def run(table_name, expected_row_count, min_row_count, max_row_count):
    sql = 'SELECT COUNT(*) FROM {}'.format(table_name)

    with engine.connect() as conn:
        results = conn.execute(text(sql))
        row = results.fetchone()
        count = row[0]

    if expected_row_count is not None:
        if count != int(expected_row_count):
            raise ValueError('Unexpected row count {}. Expected {}.'.format(str(count), str(expected_row_count)))

    if min_row_count is not None:
        if count < int(min_row_count):
            raise ValueError('Unexpected row count {}. Expected more than {}.'.format(str(count), str(min_row_count)))

    if max_row_count is not None:
        if count > int(max_row_count):
            raise ValueError('Unexpected row count {}. Expected less than {}.'.format(str(count), str(max_row_count)))

    print('Observed expected row count of {}'.format(expected_row_count))


def main():

    # Example:
    #   python -m moldb.count_rows --table-name foo --expected-rows 100

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Filter')
    parser.add_argument('-t', '--table-name', help="DB table name")
    parser.add_argument('-r', '--expected-rows', type=int, help="Expected row count")
    parser.add_argument('--min-rows', type=int, help="Minimum row count")
    parser.add_argument('--max-rows', type=int, help="Maximum row count")


    args = parser.parse_args()
    DmLog.emit_event("count_rows: ", args)

    run(args.table_name, args.expected_rows, args.min_rows, args.max_rows)


if __name__ == "__main__":
    main()
