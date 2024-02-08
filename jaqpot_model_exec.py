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


import argparse, os, sys
from jaqpotpy.models import MolecularModel
import utils
import rdkit_utils
from dm_job_utilities.dm_log import DmLog
from jaqpotpy import Jaqpot


def execute_jaqpot_model(input, output, model_id, api_key, filter=False, threshold=None, descending=False, delimiter=None,
                         id_column=None, read_header=False, write_header=False, sdf_read_records=100, interval=0):

    from jaqpotpy.cfg import config
    config.verbose = False

    jaqpot = Jaqpot()
    jaqpot.set_api_key(api_key)
    molmod = MolecularModel.load_from_jaqpot(jaqpot=jaqpot, id=model_id)
    if molmod is None:
        msg = "Model not found " + model_id
        DmLog.emit_event(msg)
        raise ValueError(msg)

    # molmod.Y can be an array in multi prediction
    model_name = molmod.Y.replace(' ', '_').replace('__', '_')
    model_title = molmod.model_title
    model_type = molmod.modeling_task
    DmLog.emit_event("Model title:", model_title, ", property name", model_name, ", type:", model_type)

    reader = rdkit_utils.create_reader(input, delimiter=delimiter, read_header=read_header,
                                       id_column=id_column, sdf_read_records=sdf_read_records)
    extra_field_names = reader.get_extra_field_names()

    calc_prop_names = get_calc_prop_names(molmod, model_name)

    if output:
        utils.expand_path(output)
        writer = rdkit_utils.create_writer(output, extra_field_names=extra_field_names, calc_prop_names=calc_prop_names,
                                           delimiter=delimiter)

    num_outputs = 0
    count = 0
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break
        count += 1
        mol, smi, id, props = t

        if count == 1 and write_header:
            headers = rdkit_utils.generate_header_values(
                reader.get_mol_field_name(), reader.field_names, len(props), calc_prop_names)
            if output:
                writer.write_header(headers)
            else:
                print(" ".join(headers))

        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records".format(count))
        if count % 10000 == 0:
            # Emit a 'total' cost, replacing all prior costs
            DmLog.emit_cost(count)

        molmod(mol)

        if filter and filter_molecule(molmod, threshold, descending):
            continue

        num_outputs += 1

        values = get_calc_values(molmod)
        if output:
            writer.write(smi, mol, id, props, values)
        else:
            print(smi, id, *values)

    DmLog.emit_event(num_outputs, "outputs among", count, "molecules")
    DmLog.emit_cost(count)


def get_calc_prop_names(molmod, prefix):
    """
    Get the names of the properties that will be output.
    These will be used as the field names of the values obtained from the get_calc_values method.
    :param molmod: The Jaqpot model to get the data from
    :param prefix: The prefix for the field names
    :return: List of names
    """
    model_type = molmod.modeling_task
    names = []
    if model_type == 'classification':
        names.append(prefix + '_Prediction')  # 0 or 1
        names.append(prefix + '_Inactive')    # probability of being inactive
        names.append(prefix + '_Active')      # probability of being active
        if molmod.doa:
            names.append(prefix + '_DOA')     # True or False
    elif model_type == 'regression':
        names.append(prefix + '_Prediction')  # regression score
        if molmod.doa:
            names.append(prefix + '_DOA')     # True or False

    return names


def get_calc_values(molmod):
    """
    Get the values that should be output.
    This depends on the type of the model.
    :param molmod: The Jaqpot model
    :return: List of values
    """
    model_type = molmod.modeling_task
    values = []
    if model_type == 'classification':
        values.append(molmod.prediction[0])
        values.append(molmod.probability[0][0])
        values.append(molmod.probability[0][1])
        if molmod.doa:
            values.append(molmod.doa.IN[0])
    elif model_type == 'regression':
        values.append(molmod.prediction[0][0])
        if molmod.doa:
            values.append(molmod.doa.IN[0])

    return values


def filter_molecule(molmod, threshold, descending):
    """
    Whether to filter this record out
    :param molmod: The Jaqpot model
    :param threshold: Filter threshold
    :param descending: Smaller regression numbers are better
    :return: If True then the molecule should be excluded
    """
    model_type = molmod.modeling_task
    if model_type == 'classification':
        if threshold is not None:
            if threshold >= molmod.probability[0][1]:
                return True
        elif not molmod.prediction[0]:
            return True
    elif model_type == 'regression' and threshold is not None:
        if threshold is not None:
            if descending:
                if threshold <= molmod.prediction[0][0]:
                    return True
            else:
                if threshold >= molmod.prediction[0][0]:
                    return True

    return False


def main():

    # command line args definitions #########################################
    parser = argparse.ArgumentParser(description='Jaqpot Model Execution')
    parser.add_argument('-m', '--model-id', required=True, help="Jaqpot model ID") # e.g. qbqUnF08SU1MhVnj3Bwd
    parser.add_argument('-i', '--input', required=True, help="Molecules to predict (.sdf or .smi)")
    parser.add_argument('-o', '--output', help="Output file (.sdf or .smi)")
    parser.add_argument('-k', "--api-key", help="Jaqpot API key")
    parser.add_argument('-f', "--filter", action='store_true', help="Only include actives in the output")
    parser.add_argument('-t', "--threshold", type=float,
                        help="Filtering threshold. If not specified then active category is used.")
    parser.add_argument('-a', "--descending", action='store_true',
                        help="Lower regression scores are better. If not set then assumed to be ascending.")
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
    DmLog.emit_event("Jaqpot model exec Args: ", args)

    delimiter = utils.read_delimiter(args.delimiter)

    api_key = None
    if args.api_key:
        api_key = args.api_key
    else:
        api_key = os.getenv('JAQPOT_API_KEY')
    if not api_key:
        DmLog.emit_event('WARNING: no Jaqpot API key provided')
        sys.exit(1)

    print('API KEY:', api_key)

    execute_jaqpot_model(args.input, args.output, args.model_id, api_key,
                         filter=args.filter, threshold=args.threshold, descending=args.descending,
                         delimiter=delimiter, id_column=args.id_column,
                         read_header=args.read_header, write_header=args.write_header,
                         sdf_read_records=args.sdf_read_records, interval=args.interval)


if __name__ == "__main__":
    main()
