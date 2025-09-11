import argparse
import glob
import shutil

from dm_job_utilities.dm_log import DmLog


def find_files(files_glob):
    files = glob.glob(files_glob)
    DmLog.emit_event("Found {} files using {}".format(len(files), files_glob))
    return files


def concat_binary(files_glob, output):
    files = find_files(files_glob)
    with (open(output, 'wb') as outfile):
        file_count = 0
        for file in files:
            file_count += 1
            with open(file,'rb') as infile:
                shutil.copyfileobj(infile, outfile)

        DmLog.emit_event("Wrote {} files".format(file_count))


def concat_text(files_glob, header, output):
    files = find_files(files_glob)
    output_count = 0
    with (open(output, 'w') as outfile):
        file_count = 0
        for file in files:
            file_count += 1

            with open(file) as infile:
                line_count = 0
                for line in infile:
                    line_count += 1
                    if header is None \
                            or (header == 'ignore' and line_count > 1) \
                            or (header == 'retain' and line_count == 1 and file_count == 1) \
                            or (header == 'retain' and line_count > 1):
                        outfile.write(line)
                        output_count += 1

    DmLog.emit_event("Wrote {} lines from {} files".format(output_count, file_count))


def main():

    # Examples:
    #   python -m concatenator -f "*.sdf"
    #   python -m concatenator -f "abcd*/output.sdf"
    #   python -m concatenator -f "*.smi" --header ignore
    #   python -m concatenator -f "*.bin" --binary
    #
    # NOTE: that if using globs for the files argument this must be escaped (e.g. abcd\*) or put in
    # quotes (e.g. "abcd*") so that they are not expanded by the shell.
    # NOTE: when using the --binary argument the --header argument is ignored.

    # command line args definitions #########################################
    parser = argparse.ArgumentParser(description='Concatenate files')
    parser.add_argument('-f', '--files', required=True, help="Name(s) of files to look for (glob allowed)")
    parser.add_argument('-o', '--output', required=True, help="Name(s) of output file")
    parser.add_argument('--header', choices=["ignore", "retain"],
                        help="Files have a header line, and what to do with it. If 'retain' the header of the first file is retained")
    parser.add_argument('-b', '--binary', action='store_true', help='Treat files as having binary content')

    args = parser.parse_args()
    DmLog.emit_event("Concatenate files: ", args)

    if args.binary:
        concat_binary(args.files, args.output)
    else:
        concat_text(args.files, args.header, args.output)


if __name__ == "__main__":
    main()
