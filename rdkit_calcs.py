#!/usr/bin/env python

# Copyright 2024 Informatics Matters Ltd.
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

import gzip, math, os, pickle
from sigfig import round

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen

import rdkit_utils


def calc_hac(mol, name=None):
    val = mol.GetNumHeavyAtoms()
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val


def calc_num_rot_bonds(mol, name=None):
    val = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val


def calc_num_rings(mol, name=None):
    val = rdMolDescriptors.CalcNumRings(mol)
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val


def calc_num_aro_rings(mol, name=None):
    val = rdMolDescriptors.CalcNumAromaticRings(mol)
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val


def calc_num_cc(mol,  name_total=None, name_undef=None):
    num_cc, num_undef_cc = rdkit_utils.get_num_chiral_centers(mol)
    if name_total is not None and num_cc is not None:
        mol.SetIntProp(name_total, num_cc)
    if name_undef is not None and num_undef_cc is not None:
        mol.SetIntProp(name_undef, num_undef_cc)
    return num_cc, num_undef_cc


def calc_num_sp3(mol, name=None):
    val = rdkit_utils.get_num_sp3_centres(mol)
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val


def calc_logp_crippen(mol, name=None, sigfigs=3):
    logp = Crippen.MolLogP(mol)
    # logp values have a silly number of decimal places
    if logp is not None:
        val = round(logp, sigfigs=sigfigs)
        if name is not None:
            mol.SetDoubleProp(name, val)
    else:
        val = None
    return val


def calc_tpsa(mol, name=None, sigfigs=3):
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    # tpsa values have a silly number of decimal places
    if tpsa is not None:
        val = round(tpsa, sigfigs=sigfigs)
        if name is not None:
            mol.SetDoubleProp(name, val)
    else:
        val = None
    return val


def calc_num_spiro(mol, name=None):
    val = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val


def calc_num_bridgehead_atoms(mol, name=None):
    val = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    if name is not None and val is not None:
        mol.SetIntProp(name, val)
    return val

########### start of sa score stuff
_sa_fscores = None


def _read_sa_fragment_scores(name='fpscores'):

    global _sa_fscores

    # generate the full path filename:
    if name == "fpscores":
        name = os.path.join(os.path.dirname(__file__), name)
    data = pickle.load(gzip.open('%s.pkl.gz' % name))
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _sa_fscores = outDict


def calc_sa_score(m, name=None):
    """
    Synthetic accessibility score.
    This is based on the work of Peter Ertl and Greg Landrum that can be found here:
    https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score

    That in turn is based on this paper:
    Peter Ertl and Ansgar Schuffenhauer
    Journal of Cheminformatics 1:8 (2009)
    http://www.jcheminf.com/content/1/1/8
    """
    if _sa_fscores is None:
        _read_sa_fragment_scores()

    # fragment score
    fp = rdMolDescriptors.GetMorganFingerprint(m, 2)  # <- 2 is the *radius* of the circular fingerprint
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0
    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += _sa_fscores.get(sfp, -4) * v
    score1 /= nf

    # features score
    nAtoms = m.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
    ri = m.GetRingInfo()
    nSpiro = calc_num_spiro(m)
    nBridgeheads = calc_num_bridgehead_atoms(m)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = 0.
    # ---------------------------------------
    # This differs from the paper, which defines:
    #  macrocyclePenalty = math.log10(nMacrocycles+1)
    # This form generates better results when 2 or more macrocycles are present
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    # correction for the fingerprint density
    # not in the original publication, added in version 1.1
    # to make highly symmetrical molecules easier to synthetise
    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.
    # smooth the 10-end
    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    if name is not None:
        m.SetDoubleProp(name, sascore)
    return sascore

########### end of sa score stuff