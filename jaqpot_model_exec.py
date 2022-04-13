import argparse, os, sys
from jaqpotpy.models import MolecularModel
import utils
import rdkit_utils
from dm_job_utilities.dm_log import DmLog
from jaqpotpy import Jaqpot


def execute_jaqpot_model(input, output, model_id, api_key, filter=False, delimiter=None, id_column=None,
                         read_header=False, write_header=False, sdf_read_records=100, interval=0):

    from jaqpotpy.cfg import config
    config.verbose = False

    jaqpot = Jaqpot()
    jaqpot.set_api_key(api_key)
    molmod = MolecularModel.load_from_jaqpot(jaqpot=jaqpot, id=model_id)
    model_name = molmod.Y.replace(' ', '_')
    DmLog.emit_event("Model is", model_name)

    reader = rdkit_utils.create_reader(input, delimiter=delimiter, read_header=read_header, sdf_read_records=sdf_read_records)
    extra_field_names = reader.get_extra_field_names()

    calc_prop_names = [model_name + '_Inactive', model_name + '_Active']

    utils.expand_path(output)
    writer = rdkit_utils.create_writer(output, extra_field_names=extra_field_names, calc_prop_names=calc_prop_names,
                                       delimiter=delimiter)
    if write_header:
        headers = ['smiles']
        if extra_field_names:
            headers.extend(extra_field_names)
        else:
            for i, prop in enumerate(calc_prop_names):
                headers.append('field' + str(i + 2))
        headers.extend(calc_prop_names)
        writer.write_header(headers)

    num_actives = 0
    count = 0
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break
        count += 1
        mol, smi, id, props = t

        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records".format(count))
        if count % 10000 == 0:
            # Emit a 'total' cost, replacing all prior costs
            DmLog.emit_cost(count)

        molmod(mol)
        if filter and not molmod.prediction[0]:
            continue
        utils.log("mol:", count, molmod.prediction, molmod.probability)
        writer.write(smi, mol, id, props, (molmod.probability[0]))
        if molmod.prediction[0]:
            num_actives += 1

    DmLog.emit_event("Found", num_actives, "actives among", count, "molecules")
    DmLog.emit_cost(count)


def main():

    # command line args definitions #########################################
    parser = argparse.ArgumentParser(description='Jaqpot Model Execution')
    parser.add_argument('-m', '--model-id', required=True, help="Jaqpot model ID") # e.g. qbqUnF08SU1MhVnj3Bwd
    parser.add_argument('-i', '--input', required=True, help="Molecules to predict (.sdf or .smi)")
    parser.add_argument('-o', '--output', required=True, help="Output file (.sdf or .smi)")
    parser.add_argument('-k', "--api-key", help="Jaqpot API key")
    parser.add_argument('-f', "--filter", action='store_true', help="Only include actives in the output")
    # to pass tab as delimiter specify it as $'\t' or use one of the symbolic names 'comma', 'tab', 'space' or 'pipe'
    parser.add_argument('-d', '--delimiter', help="Delimiter when using SMILES")
    parser.add_argument('--id-column', help="Column for name field (zero based integer for .smi, text for SDF)")
    parser.add_argument('--read-header', action='store_true',
                        help="Read a header line with the field names when reading .smi or .txt")
    parser.add_argument('--write-header', action='store_true', help='Write a header line when writing .smi or .txt')
    parser.add_argument('--sdf-read-records', default=100, type=int,
                        help="Read this many SDF records to determine field names")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("Expand Args: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    api_key = None
    if args.api_key:
        api_key = args.api_key
    else:
        api_key = os.getenv('JAQPOT_API_KEY')
    if not api_key:
        DmLog.emit_event('WARNING: no Jaqpot API key provided')
        sys.exit(1)

    execute_jaqpot_model(args.input, args.output, args.model_id, api_key, filter=args.filter, delimiter=delimiter,
                         id_column=args.id_column,
                         read_header=args.read_header, write_header=args.write_header,
                         sdf_read_records=args.sdf_read_records, interval=args.interval)


if __name__ == "__main__":
    main()
