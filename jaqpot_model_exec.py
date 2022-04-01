import argparse
from jaqpotpy.models import MolecularModel
from rdkit import Chem
import utils
import rdkit_utils
from jaqpotpy import Jaqpot


def execute(input, output, model_id, api_key):

    from jaqpotpy.cfg import config
    config.verbose = False

    reader = rdkit_utils.create_reader(input)
    extra_field_names = reader.get_extra_field_names()

    utils.expand_path(output)
    # writer = rdkit_utils.create_writer(output, extra_field_names=extra_field_names)

    mols = []

    jaqpot = Jaqpot()
    jaqpot.set_api_key(api_key)
    molmod = MolecularModel.load_from_jaqpot(jaqpot=jaqpot, id=model_id)
    suppl = Chem.SDMolSupplier(input)
    for mol in suppl:
        molmod(mol)
        print(molmod.Y)
        print(molmod.prediction)
        print(molmod.doa)
        print(molmod.probability)


def main():

    # command line args definitions #########################################
    parser = argparse.ArgumentParser(description='Jaqpot Model')
    parser.add_argument('-m', '--model-id', required=True, help="Jaqpot model ID") # e.g. qbqUnF08SU1MhVnj3Bwd
    parser.add_argument('-i', '--input', required=True, help="Molecules to predict (.sdf or .smi)")
    parser.add_argument('-o', '--output', required=True, help="Output file (.sdf or .smi)")
    parser.add_argument('-k', "--api-key", help="Jaqpot API key")

    args = parser.parse_args()

    api_key = None
    if args.api_key:
        api_key = args.api_key
    else:
        api_key = os.getenv('JAQPOT_API_KEY')
    if not api_key:
        utils.log('WARNING: no Jaqpot API key provided')

    execute(args.input, args.output, args.model_id, api_key)



if __name__ == "__main__":
    main()
