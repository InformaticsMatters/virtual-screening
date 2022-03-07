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
Basic SuCOS scoring. Allows a set of molecules from a SD file to be overlayed to a reference molecule,
with the resulting scores being written as properties in the output SD file.

SuCOS is the work of Susan Leung.
GitHub: https://github.com/susanhleung/SuCOS
Publication: https://doi.org/10.26434/chemrxiv.8100203.v1
"""


import argparse, os
import numpy as np
from rdkit import rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
import utils
import rdkit_utils
from dm_job_utilities.dm_log import DmLog


### start function definitions #########################################

field_SuCOS_Score = "SuCOS_Score"
field_SuCOS_FMScore = "SuCOS_FeatureMap_Score"
field_SuCOS_TaniScore = "SuCOS_Tanimoto_Score"
field_SuCOS_ProtrudeScore = "SuCOS_Protrude_Score"

# Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')


def filter_feature(f):
    result = f.GetFamily() in keep
    # TODO - nothing ever seems to be filtered. Is this expected?
    if not result:
        utils.log("Filtered out feature type", f.GetFamily())
    return result


def get_raw_features(mol):
    raw_feats = fdef.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're interested in
    filtered = list(filter(filter_feature, raw_feats))
    return filtered


def get_feature_map_score(small_feats, large_feats, tani=False, score_mode=FeatMaps.FeatMapScoreMode.All):
    """
    Generate the feature map score.

    :param small_feats:
    :param large_feats:
    :param tani:
    :return:
    """

    featLists = []
    for rawFeats in [small_feats, large_feats]:
        # filter that list down to only include the ones we're interested in
        featLists.append(rawFeats)
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    # set the score mode
    fms[0].scoreMode = score_mode

    try:
        if tani:
            c = fms[0].ScoreFeats(featLists[1])
            A = fms[0].GetNumFeatures()
            B = len(featLists[1])
            if B != fms[1].GetNumFeatures():
                utils.log("Why isn't B equal to number of features...?!")
            tani_score = float(c) / (A+B-c)
            return tani_score
        else:
            fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
            return fm_score
    except ZeroDivisionError:
        utils.log("ZeroDivisionError")
        return 0.0

    if tani:
        tani_score = float(c) / (A+B-c)
        return tani_score
    else:
        fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
        return fm_score


def get_sucos_scores(ref_mol, query_mol, tani=False, ref_features=None, query_features=None, score_mode=FeatMaps.FeatMapScoreMode.All):
    """
    This is the key function that calculates the SuCOS scores and is expected to be called from other modules.
    To improve performance you can pre-calculate the features and pass them in as optional parameters to avoid having
    to recalculate them. Use the get_raw_features function to pre-calculate the features.

    :param ref_mol: The reference molecule to compare to
    :param query_mol: The molecule to compare to the reference
    :param tani: Whether to calculate Tanimoto distances
    :param ref_features: An optional feature map for the reference molecule, avoiding the need to re-calculate it.
    :param query_features: An optional feature map for the query molecule, avoiding the need to re-calculate it.
    :return: A tuple of 3 values. 1 the sucos score, 2 the feature map score,
        3 the Tanimoto distance or 1 minus the protrude distance
    """

    if not ref_features:
        ref_features = get_raw_features(ref_mol)
    if not query_features:
        query_features = get_raw_features(query_mol)

    try:
        fm_score = get_feature_map_score(ref_features, query_features, tani, score_mode)
        fm_score = np.clip(fm_score, 0, 1)

        if tani:
            tani_sim = 1 - float(rdShapeHelpers.ShapeTanimotoDist(ref_mol, query_mol))
            tani_sim = np.clip(tani_sim, 0, 1)
            sucos_score = 0.5*fm_score + 0.5*tani_sim
            return (sucos_score, fm_score, tani_sim)
        else:
            protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref_mol, query_mol, allowReordering=False)
            protrude_dist = np.clip(protrude_dist, 0, 1)
            protrude_val = 1.0 - protrude_dist
            sucos_score = 0.5 * fm_score + 0.5 * protrude_val
            return (sucos_score, fm_score, protrude_val)

    except Exception:
        utils.log("Failed to calculate SuCOS scores. Returning 0,0,0")
        return 0.0, 0.0, 0.0


def process(reference, input, output, tani=False, score_mode=FeatMaps.FeatMapScoreMode.All, interval=None):

    ref_mol = rdkit_utils.rdk_read_single_mol(reference)
    DmLog.emit_event("Target mol has {} heavy atoms".format(ref_mol.GetNumHeavyAtoms()))
    ref_features = get_raw_features(ref_mol)

    # create reader
    calc_prop_names = [field_SuCOS_Score, field_SuCOS_FMScore]
    if tani:
        calc_prop_names.append(field_SuCOS_TaniScore)
    else:
        calc_prop_names.append(field_SuCOS_ProtrudeScore)

    reader = rdkit_utils.SdfReader(input, 0, 50)

    # create writer
    utils.expand_path(output)
    writer = rdkit_utils.SdfWriter(output, calc_prop_names)

    count = 0
    errors = 0
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break

        mol, smi, id, props = t
        count += 1
        if interval and count % interval == 0:
            DmLog.emit_event("Processed {} records".format(count))

        if not mol:
            utils.log('Failed to read molecule', count)
            errors += 1
            continue
        try:
            scores = get_sucos_scores(ref_mol, mol, tani=tani, ref_features=ref_features, score_mode=score_mode)
            # utils.log("Scores:", scores[0], scores[1], scores[2])
            writer.write(None, mol, id, props, scores)
        except ValueError as e:
            errors += 1
            utils.log("Molecule", count, "failed to score:", e.message)

    writer.close()
    DmLog.emit_event("Completed. Processed {} molecules, {} errors".format(count, errors))

    return count, errors


def parse_score_mode(value):
    if value is None or value == 'all':
        return FeatMaps.FeatMapScoreMode.All
    elif value == 'closest':
        return FeatMaps.FeatMapScoreMode.Closest
    elif value == 'best':
        return FeatMaps.FeatMapScoreMode.Best
    else:
        raise ValueError(value + " is not a valid scoring mode option")


### start main execution #########################################

def main():

    parser = argparse.ArgumentParser(description='SuCOS with RDKit')

    parser.add_argument('-i', '--input', required=True, help="File with molecules to cluster (.sdf)")
    parser.add_argument('-o', '--output', required=True, help="Output file (.sdf)")
    parser.add_argument('-r', '--reference', help='Target molecule to compare against (.sdf or .mol')
    parser.add_argument('-t', '--tanimoto', action='store_true', help='Include Tanimoto distance in score')
    parser.add_argument('-m', '--score-mode', choices=['all', 'closest', 'best'],
                        help="choose the scoring mode for the feature map, default is 'all'.")
    parser.add_argument("--interval", type=int, help="Reporting interval")

    args = parser.parse_args()
    DmLog.emit_event("SuCOS Args: ", args)

    score_mode = parse_score_mode(args.score_mode)

    # this does the processing
    count, errors = process(args.reference, args.input, args.output,
                            tani=args.tanimoto, score_mode=score_mode, interval=args.interval)


if __name__ == "__main__":
    main()
