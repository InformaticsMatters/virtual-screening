import argparse
import glob
import shutil

from dm_job_utilities.dm_log import DmLog


def find_files(input_file, dirs_glob):

    files = glob.glob(f"{dirs_glob}/{input_file}")
    DmLog.emit_event("Found {} files using {}".format(len(files), files_glob))
    return files


def concat_binary(input_file, dirs_glob, output):
    files = find_files(input_file, dirs_glob)
    with (open(output, 'wb') as outfile):
        file_count = 0
        for file in files:
            file_count += 1
            with open(file,'rb') as infile:
                shutil.copyfileobj(infile, outfile)

        DmLog.emit_event("Wrote {} files".format(file_count))


def concat_text(input_file, dirs_glob, header, output):
    files = find_files(input_file, dirs_glob)
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
    #   python -m concatenator -f "output.sdf" -d "input-*"
    #
    # NOTE: when using globs for the files argument this must be escaped (e.g. abcd\*) or put in
    # quotes (e.g. "abcd*") so that they are not expanded by the shell.
    # NOTE: when using the --binary argument the --header argument is ignored.

    # command line args definitions #########################################
    parser = argparse.ArgumentParser(description='Concatenate files')
    parser.add_argument('-f', '--input-file', required=True, help="Name of the file to concatenate")
    parser.add_argument('-d', '--dirs-glob', required=True, help="Glob of directories to search")
    parser.add_argument('-o', '--output', required=True, help="Name(s) of output file")
    parser.add_argument('--header', choices=["ignore", "retain"],
                        help="Files have a header line, and what to do with it. If 'retain' the header of the first file is retained")
    parser.add_argument('-b', '--binary', action='store_true', help='Treat files as having binary content')

    args = parser.parse_args()
    DmLog.emit_event("Concatenate files: ", args)

    if args.binary:
        concat_binary(args.input_file, args.dirs_glob, args.output)
    else:
        concat_text(args.input_file, args.dirs_glob, args.header, args.output)


if __name__ == "__main__":
    main()
