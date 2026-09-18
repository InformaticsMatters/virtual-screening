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

"""
Job that allows failures to be tested
"""

import argparse
from dm_job_utilities.dm_log import DmLog

def main():
    # Example usage:
    # ./fail.py --mode error|empty|missing

    parser = argparse.ArgumentParser(description='Fail')
    parser.add_argument('-m', '--mode', required=True, choices=['error', 'empty', 'missing'], help='Fail mode')
    parser.add_argument('-o', '--output', help='Output file name')

    args = parser.parse_args()
    DmLog.emit_event("fail: ", args)

    if args.output:
        outfile = args.output
    else:
        outfile = 'test-fail-results.txt'

    if args.mode == 'error':
        DmLog.emit_event('Job failed')
        exit(1)
    elif args.mode == 'empty':
        with open(outfile, 'wt') as out:
            DmLog.emit_event('Empty test-fail-results.txt')
    elif args.mode == 'missing':
        DmLog.emit_event('Missing test-fail-results.txt')

if __name__ == "__main__":
    main()
