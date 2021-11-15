#!/usr/bin/env python
"""Convert from one file type to another file based on the mime-type

    Supported conversions:
    - sdf to json.

"""
import argparse
import json
import uuid
import os
import gzip
import sys
from typing import Dict
from conversion_utils import sdf_get_next_record, is_valid_uuid
from utils import log_dm_event

_CONVERSION_MAP = {'chemical/x-mdl-sdfile': 'sdf',
                  'squonk/x-dataset-molecule-v2+json': 'json'}

# Main Class
class ConvertFile:
    """Class ConvertFile

    Purpose: Converts from an input mime-type to an output mime-type.

    """
    errors: int = 0
    lines: int = 0
    records: int = 0
    interval: int = 0

    def __init__(self, interval: int = 0):
        self.interval = interval

    def convert(self,
                from_mime_type: str,
                to_mime_type: str,
                infile: str,
                outfile: str) -> bool:

        """Dispatch method"""
        from_type: str = _CONVERSION_MAP.get(from_mime_type)
        to_type: str = _CONVERSION_MAP.get(to_mime_type)

        method_name = 'convert_' + str(from_type) + '_to_' + str(to_type)
        # Get the method from 'self'.
        try:
            method = getattr(self, method_name)
        except AttributeError:
            log_dm_event('Method to support %s not found', to_mime_type)
            return False

        # Call the method as we return it
        return method(infile, outfile) # pylint: disable=too-many-function-args

    def process_molecules_json(self, infile, outfile):
        """ process molecules in SDF file
        """

        # Loop through molecules
        while True:
            molecule: Dict = {}
            molecule_block, molecule_name, properties =\
                sdf_get_next_record(infile)

            if not molecule_block:
                break

            if self.records:
                outfile.write(',')
            self.records += 1

            # title line
            if molecule_name:
                molecule['name'] = molecule_name

            if is_valid_uuid(molecule_name):
                molecule_uuid = molecule_name
            else:
                molecule_uuid = str(uuid.uuid4())

            molecule['molblock'] = molecule_block

            record = {'uuid': molecule_uuid, 'molecule': molecule,
                      'values': properties}
            json_str = json.dumps(record)
            outfile.write(json_str)

            if self.interval and self.records % self.interval == 0:
                log_dm_event("Processed {} records".format(self.records))

        outfile.write(']')

    def convert_sdf_to_json(self, infile, outfile) -> bool:
        """converts the given SDF file into a Squonk json file.
           Returns True if file successfully converted.
        """

        if infile.endswith('.gz'):
            infile_handle = gzip.open(infile, 'rt')
        else:
            infile_handle = open(infile, 'rt')

        if outfile.endswith('.gz'):
            outfile_handle = gzip.open(outfile, 'wt')
        else:
            outfile_handle = open(outfile, 'wt')

        outfile_handle.write('[')

        try:
            self.process_molecules_json(infile_handle, outfile_handle)
        finally:
            outfile_handle.close()

        if self.errors > 0:
            return False
        return True


def main():
    # Example usage:
    #   python convert_file.py -i data/test.sdf -if chemical/x-mdl-sdfile\
    #     -o test.json -of squonk/x-dataset-molecule-v2+json --interval 10000

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Convert file')
    parser.add_argument('-i', '--input-file', required=True,
                        help="Input file")
    parser.add_argument('-if', '--input-format', required=True,
                        help="Input file format (mime-type)")
    parser.add_argument('-o', '--output-file', required=True,
                        help="Output file")
    parser.add_argument('-of', '--output-format', required=True,
                        help="Output file format (mime-type")
    parser.add_argument("--interval", type=int,
                        help="Reporting interval")

    args = parser.parse_args()
    log_dm_event("convert_file.py: ", args)

     # if no input file then raise error and exit
    if not os.path.isfile(args.input_file):
        log_dm_event('File {} is not present'.format(args.input_file))
        sys.exit(1)

    log_dm_event('Converting {} to format {} in file {}...'.
                 format(args.input_file, args.output_format, args.output_file))

    converter = ConvertFile(args.interval)

    processed: bool = converter.convert(args.input_format,
                                        args.output_format,
                                        args.input_file,
                                        args.output_file)

    if processed:
        log_dm_event('Converter finished successfully')
        log_dm_event('records processed={}'.format(converter.records))
        log_dm_event('errors={}'.format(converter.errors))
    else:
        log_dm_event('Converter failed')
        sys.exit(1)


if __name__ == "__main__":
    main()
