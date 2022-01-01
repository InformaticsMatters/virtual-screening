#!/usr/bin/env python

# Copyright 2020 Informatics Matters Ltd.
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
Classes and utility functions for calculating interactions using ODDT.
"""

import sys, json, math

import numpy as np

import oddt
from oddt import toolkit, spatial, interactions

I_TYPE_HBOND = 'HydrogenBond'
I_TYPE_HALOGEN = 'HalogenBond'
I_TYPE_HYDROPHOBIC = 'Hydrophobic'
I_TYPE_SALT_BRIDGE = 'SaltBridge'
I_TYPE_PI_STACKING = 'PiStacking'
I_TYPE_PI_CATION = 'PiCation'

I_SUFFIX = 'Interaction'
I_NAME_HBOND = I_TYPE_HBOND + I_SUFFIX
I_NAME_HALOGEN = I_TYPE_HALOGEN + I_SUFFIX
I_NAME_HYDROPHOBIC = I_TYPE_HYDROPHOBIC + I_SUFFIX
I_NAME_SALT_BRIDGE = I_TYPE_SALT_BRIDGE + I_SUFFIX
I_NAME_PI_STACKING = I_TYPE_PI_STACKING + I_SUFFIX
I_NAME_PI_CATION = I_TYPE_PI_CATION + I_SUFFIX


class Interaction:
    def __init__(self, canon_residue, protein_pos, ligand_pos, distance, ligand_atom):
        self.canon_residue = canon_residue
        self.protein_pos = protein_pos
        self.ligand_pos = ligand_pos
        self.distance = distance
        self.ligand_atom = ligand_atom
        self.score_value = None
        self.score_name = None

    def compare(self, other):
        p_dist = distance_between_points(self.protein_pos, other.protein_pos)
        l_dist = distance_between_points(self.ligand_pos, other.ligand_pos)
        score = 0
        # score += max(0.0, (2.0 - p_dist) / 4.0)
        score += max(0.0, (2.0 - l_dist) / 2.0)
        # if score > 0:
        #     print('       ', score, l_dist, self.ligand_pos, other.ligand_pos)

        return score

    def compare_and_store(self, other, name):
        sc = self.compare(other)
        if not self.score_value or sc > self.score_value:
            self.score_value = sc
            self.score_name = name
        return sc

    def __eq__(self, other):
        if not isinstance(other, Interaction):
            return False
        else:
            return self.canon_residue == other.canon_residue and self.protein_pos == other.protein_pos and self.ligand_pos == other.ligand_pos

    def __hash__(self):
        return hash(self.canon_residue, self.protein_pos, self.ligand_pos)


class InteractionType:
    def __init__(self, interaction_type, interactions):
        self.interaction_type = interaction_type
        self.interactions = interactions

    def addInteraction(self, interaction):
        self.interactions.append(interaction)

    def compare(self, other):
        if self.interaction_type != other.interaction_type:
            print('WARNING: comparing interactions of different types:', self.interaction_type, other.interaction_type)
        scores = []
        data = []
        for inter1 in self.interactions:
            for inter2 in other.interactions:
                if inter1.canon_residue == inter2.canon_residue:
                    score = inter1.compare(inter2)
                    scores.append(score)
                    if score > 0:
                        # print('  ', self.interaction_type, inter1.canon_residue, score)
                        data.append((inter1.canon_residue, score))
        return sum(scores), len(scores), data

    def asText(self):
        content = []
        for i in self.interactions:
            if i.ligand_atom is not None:
                s = "%s, [%s, %s, %s], [%s, %s, %s], [%.3f, %s]" % (
                i.canon_residue, i.protein_pos[0], i.protein_pos[1], i.protein_pos[2],
                i.ligand_pos[0], i.ligand_pos[1], i.ligand_pos[2],
                i.distance, i.ligand_atom)
            else:
                s = "%s, [%s, %s, %s], [%s, %s, %s], [%.3f]" % (
                i.canon_residue, i.protein_pos[0], i.protein_pos[1], i.protein_pos[2],
                i.ligand_pos[0], i.ligand_pos[1], i.ligand_pos[2], i.distance)
            if i.score_value:
                s += ", [%s, %s]" % (i.score_value, i.score_name)
            content.append(s)
        return "\n".join(content)


class InteractionSet:
    def __init__(self, id, ligand_name, interaction_types):
        self.id = id
        self.ligand_name = ligand_name
        self.interaction_types = interaction_types

    def add(self, interaction_type):
        self.interaction_types.append(interaction_type)

    def compare(self, other):
        scores = []
        matches = []
        for itype1 in self.interaction_types:
            for itype2 in other.interaction_types:
                if itype1.interaction_type == itype2.interaction_type:
                    score, count, data = itype1.compare(itype2)
                    if count:
                        # print('  ', itype1.interaction_type, score, count)
                        matches.append((itype1.interaction_type, data))
                    scores.append(score)
        return sum(scores), len(scores), matches


class InteractionEncoder(json.JSONEncoder):
    def default(self, i):
        if isinstance(i, Interaction):
            return i.__dict__
        elif isinstance(i, InteractionType):
            return {i.interaction_type: i.interactions}
        elif isinstance(i, InteractionSet):
            d = {} #collections.OrderedDict()
            d['id'] = i.id
            d['ligand_name'] = i.ligand_name
            for t in i.interaction_types:
                d[t.interaction_type] = t.interactions
            return d
        else:
            return super().default(i)


def distance_between_points(pos1, pos2):
    return math.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2 + (pos1[2] - pos2[2]) ** 2)


def from_json(text):
    results = []
    records = json.loads(text)
    for record in records:
        itypes = []
        name = None
        for key in record:
            if key == 'id':
                id = record['id']
            elif key == 'ligand_name':
                name = record['ligand_name']
            elif key.endswith('Interaction'):
                itype = InteractionType(key, [])
                itypes.append(itype)
                values = record[key]
                for value in values:
                    inter = Interaction(value['canon_residue'], value['protein_pos'], value['ligand_pos'],
                                        value['distance'], value.get('ligand_atom', None))
                    itype.addInteraction(inter)
            else:
                print('WARNING: unexpected field %s' % (key))

        iset = InteractionSet(id, name, itypes)
        results.append(iset)
    return results


# def compare_interactions(reference_file, test_file):
#     with open(reference_file, "r") as ref:
#         txt = ref.read()
#         ref_data = from_json(txt)
#         print('Found', len(ref_data), 'reference items')
#     with open(test_file, "r") as test:
#         txt = test.read()
#         test_data = from_json(txt)
#         print('Found', len(test_data), 'test items')
#     for i, iset1 in enumerate(test_data, start=1):
#         canonical_sites = {}
#         for j, iset2 in enumerate(ref_data, start=1):
#             # print('Comparing', i, j, iset1.ligand_name, iset2.ligand_name)
#             score, count, matches = iset1.compare(iset2)
#             if score:
#                 for match in matches:
#                     if match[1]:
#                         # print('Matched', i, j, iset1.ligand_name, iset2.ligand_name, score)
#                         for inter in match[1]:
#                             # print('  ', match[0], inter[0], inter[1])
#                             type_canon = (match[0], inter[0])
#                             if type_canon in canonical_sites:
#                                 current_best = canonical_sites[type_canon]
#                                 if inter[1] > current_best[0][1]:
#                                     canonical_sites[type_canon] = (inter, j, iset2.ligand_name)
#                             else:
#                                 canonical_sites[type_canon] = (inter, j, iset2.ligand_name)
#         total_score = 0
#         for type_canon in canonical_sites:
#             inter_j_name = canonical_sites[type_canon]
#             total_score += inter_j_name[0][1]
#
#         print('Comparing', i, iset1.ligand_name, total_score)
#         for type_canon in canonical_sites:
#             inter_j_name = canonical_sites[type_canon]
#             print('  ', type_canon[0], type_canon[1], inter_j_name[0][1], inter_j_name[1], inter_j_name[2])


def process_mol(protein, mol, key_inters, index, filter_strict=False,
                exact_protein=False, exact_ligand=False, compare_data=None):

    mol_key_inters = []
    inter_data = {}
    num_inters = 0
    all_scores = {}

    inter_set = InteractionSet(index, None, [])

    inter_data['Index'] = index
    if '_Name' in mol.data:
        molname = mol.data['_Name']
        inter_data['LigandName'] = molname
        inter_set.ligand_name = molname

    # handle H-bond interactions
    inters = calc_hydrogen_bond_interactions(protein, mol, key_inters.get(I_TYPE_HBOND, None),
                                                      mol_key_inters, filter_strict=filter_strict,
                                                      exact_protein=exact_protein, exact_ligand=exact_ligand)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[I_NAME_HBOND] = data
        num_inters += len(inters.interactions)
        inter_data[I_NAME_HBOND] = data
        inter_set.add(inters)

    # handle hydrophobic interactions
    inters = calc_hydrophobic_interactions(protein, mol, key_inters.get(I_TYPE_HYDROPHOBIC, None),
                                                    mol_key_inters)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[I_NAME_HYDROPHOBIC] = data
        num_inters += len(inters.interactions)
        inter_data[I_NAME_HYDROPHOBIC] = data
        inter_set.add(inters)

    # handle salt bridge interactions
    inters = calc_salt_bridge_interactions(protein, mol, key_inters.get(I_TYPE_SALT_BRIDGE, None),
                                                    mol_key_inters, exact_protein=exact_protein, exact_ligand=exact_ligand)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[I_NAME_SALT_BRIDGE] = data
        num_inters += len(inters.interactions)
        inter_data[I_NAME_HYDROPHOBIC] = data
        inter_set.add(inters)

    # handle pi stacking interactions
    inters = calc_pi_stacking_interactions(protein, mol, key_inters.get(I_TYPE_PI_STACKING, None),
                                                    mol_key_inters, filter_strict=filter_strict)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[I_NAME_PI_STACKING] = data
        num_inters += len(inters.interactions)
        inter_data[I_NAME_PI_STACKING] = data
        inter_set.add(inters)

    # handle pi cation interactions
    inters = calc_pi_cation_interactions(protein, mol, key_inters.get(I_TYPE_PI_CATION, None),
                                                  mol_key_inters, filter_strict=filter_strict,
                                                  exact_protein=exact_protein, exact_ligand=exact_ligand)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[I_NAME_PI_CATION] = data
        num_inters += len(inters.interactions)
        inter_data[I_NAME_PI_CATION] = data
        inter_set.add(inters)

    # handle halogen bond interactions
    inters = calc_halogen_bond_interactions(protein, mol, key_inters.get(I_TYPE_HALOGEN, None),
                                                     mol_key_inters, filter_strict=filter_strict)
    if inters:
        compare_interactions(inters, compare_data, all_scores)
        data = inters.asText()
        mol.data[I_NAME_HALOGEN] = data
        num_inters += len(inters.interactions)
        inter_data[I_NAME_HALOGEN] = data
        inter_set.add(inters)

    mol.data['NumTotalInteractions'] = str(num_inters)
    mol.data['NumKeyInteractions'] = str(len(mol_key_inters))
    if mol_key_inters:
        keydata = []
        for key_inter in mol_key_inters:
            score = all_scores.get(key_inter, 0)
            mol.data[key_inter] = score
            keydata.append(key_inter + ' ' + str(score))
        mol.data['KeyInteractions'] = '\n'.join(keydata)

    return inter_set


def compare_interactions(inter_type, compare_to, all_scores):
    if inter_type and compare_to:
        for inter in inter_type.interactions:
            # print('  looking at', inter_type.interaction_type, inter.canon_residue, inter.ligand_pos)
            for compare_set in compare_to:
                for compare_type in compare_set.interaction_types:
                    if compare_type.interaction_type == inter_type.interaction_type:
                        # print('    analysing', compare_set.ligand_name, compare_type.interaction_type, len(compare_type.interactions))
                        for compare_inter in compare_type.interactions:
                            # print('      checking', inter.canon_residue, compare_inter.canon_residue)
                            if inter.canon_residue == compare_inter.canon_residue:
                                # print('     ', compare_inter.canon_residue, compare_set.ligand_name)
                                sc = inter.compare_and_store(compare_inter, compare_set.ligand_name)

            if inter.score_value:
                k = inter_type.interaction_type[:-11] + ':' + inter.canon_residue
                # print('Found key interaction', k, inter.score_value)
                if k in all_scores:
                    # keep the best score
                    if inter.score_value > all_scores[k]:
                        all_scores[k] = inter.score_value
                else:
                    all_scores[k] = inter.score_value


def calc_hydrogen_bond_interactions(protein, mol, key_inters_defs, mol_key_inters,
                                    filter_strict=False, exact_protein=False, exact_ligand=False):
    """ Calculate H-bond interactions

    Parameters:
    protein (Molecule): The protein
    mol (Molecule): The ligand to test
    key_inters_defs (List): Optional list of key H-bond interactions
    mol_key_inters (List): List to add the key H-bond interactions to
    filter_strict (Bool): Whether to use strict matching

    Returns:
    InteractionType: The interactions
    """

    # first apply a fix that's needed to handle the mis-assignment of donor atoms for a particular tautomeric
    # form: nc[n;H1]. See https://groups.google.com/g/oddt/c/fqzmhSqprw8/m/nmaaUlCDAgAJ
    # mol.atom_dict.setflags(write=True)
    # matches = oddt.toolkit.Smarts('[n;H0]').findall(mol)
    # for (idx,) in matches:
    #     # idx assumes 0 based indexing e.g. RDKit. OBabel uses 1 based indexing.
    #     mol.atom_dict['isdonor'][idx] = False
    # mol.atom_dict.setflags(write=False)
    # end of fix

    protein_atoms, ligand_atoms, strict = oddt.interactions.hbonds(protein, mol,
                                                                   mol1_exact=exact_protein, mol2_exact=exact_ligand)
    inters = {}

    for p, l, s in zip(protein_atoms, ligand_atoms, strict):
        if s or not filter_strict:
            c = get_canonical_hbond(p)
            dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
            p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
            l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
            t = Interaction(c, p_coords, l_coords, dist, l['id'].item())
            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(I_TYPE_HBOND + ':' + c)
    if inters:
        # print('  found', len(inters), 'h-bonds')
        return InteractionType(I_NAME_HBOND, list(inters.values()))
    else:
        return None


def calc_hydrophobic_interactions(protein, mol, key_inters_defs, mol_key_inters):
    """ Calculate hydrophobic interactions.
    ODDT generates hydrophobic interactions for all pairs of atoms that are within the specification which results in
    multiple interactions between the same hydrophobic region of the ligand and protein.
    We reduce those down to a single interaction, the one that is the shortest.
    ODDT also seems to generate hydrophobic interactions that are complete duplicates, and the de-duplication process
    removes these as well.

    Parameters:
    protein (Molecule): The protein
    mol (Molecule): The ligand to test
    key_inters_defs (List): Optional list of key H-bond interactions
    mol_key_inters (List): List to add the key H-bond interactions to

    Returns:
    InteractionType: The interactions
    """
    inters = {}
    protein_atoms, ligand_atoms = oddt.interactions.hydrophobic_contacts(protein, mol)
    for p, l in zip(protein_atoms, ligand_atoms):
        c = get_canonical_residue(p)
        dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
        p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
        l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
        t = Interaction(c, p_coords, l_coords, dist, l['id'].item())
        if c in inters:
            current_value = inters[c]
            if dist < current_value.distance:
                inters[c] = t
        else:
            inters[c] = t
            if key_inters_defs and c in key_inters_defs:
                mol_key_inters.append(I_TYPE_HYDROPHOBIC + ':' + c)

    if inters:
        # print('  found', len(inters), 'hydrophobics')
        return InteractionType(I_NAME_HYDROPHOBIC, list(inters.values()))
    else:
        return None


def calc_salt_bridge_interactions(protein, mol, key_inters_defs, mol_key_inters, exact_protein=False, exact_ligand=False):
    inters = {}
    protein_atoms, ligand_atoms = oddt.interactions.salt_bridges(protein, mol, mol1_exact=exact_protein, mol2_exact=exact_ligand)
    for p, l in zip(protein_atoms, ligand_atoms):
        c = get_canonical_residue(p)
        dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
        p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
        l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
        t = Interaction(c, p_coords, l_coords, dist, l['id'].item())
        if c in inters:
            current_value = inters[c]
            if dist < current_value.distance:
                inters[c] = t
        else:
            inters[c] = t
            if key_inters_defs and c in key_inters_defs:
                mol_key_inters.append(I_TYPE_SALT_BRIDGE + ':' + c)

    if inters:
        return InteractionType(I_NAME_SALT_BRIDGE, list(inters.values()))
    else:
        return None


def calc_pi_stacking_interactions(protein, mol, key_inters_defs, mol_key_inters, filter_strict=False):
    protein_atoms, ligand_atoms, strict_parallel, strict_perpendicular = oddt.interactions.pi_stacking(protein, mol)
    inters = {}
    for p, l, s1, s2 in zip(protein_atoms, ligand_atoms, strict_parallel, strict_perpendicular):
        if (s1 or s2) or not filter_strict:
            c = get_canonical_residue(p)
            dist = spatial.distance(np.array([p['centroid']]), np.array([l['centroid']]))[0][0]
            p_coords = (p['centroid'][0].item(), p['centroid'][1].item(), p['centroid'][2].item())
            l_coords = (l['centroid'][0].item(), l['centroid'][1].item(), l['centroid'][2].item())
            t = Interaction(c, p_coords, l_coords, dist, None)

            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(I_TYPE_PI_STACKING + ':' + c)

    if inters:
        return InteractionType(I_NAME_PI_STACKING, list(inters.values()))
    else:
        return None


def calc_pi_cation_interactions(protein, mol, key_inters_defs, mol_key_inters,
                                filter_strict=False, exact_protein=False, exact_ligand=False):
    """
    Pi-cation calculations are directional so are run in both directions (protein->ligand and ligand-> protein).
    Hence this method runs the pi_cation() calculation twice.
    :param protein:
    :param mol:
    :param key_inters_defs:
    :param mol_key_inters:
    :param filter_strict:
    :param exact_protein:
    :param exact_ligand:
    :return:
    """
    inters = {}
    # first treat the protein as the pi and the ligand as the cation
    rings, cation, strict = oddt.interactions.pi_cation(protein, mol, cation_exact=exact_ligand)
    for ring, cat, s in zip(rings, cation, strict):
        if s or not filter_strict:
            dist = spatial.distance(np.array([ring['centroid']]), np.array([cat['coords']]))[0][0]
            c = get_canonical_residue(ring)
            p_coords = (ring['centroid'][0].item(), ring['centroid'][1].item(), ring['centroid'][2].item())
            l_coords = (cat['coords'][0].item(), cat['coords'][1].item(), cat['coords'][2].item())
            t = Interaction(c, p_coords, l_coords, dist, None)

            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(I_TYPE_PI_CATION + ':' + c)

    # now with the ligand as the pi and the protein as the cation
    rings, cation, strict = oddt.interactions.pi_cation(mol, protein, cation_exact=exact_protein)
    for ring, cat, s in zip(rings, cation, strict):
        if s or not filter_strict:
            dist = spatial.distance(np.array([ring['centroid']]), np.array([cat['coords']]))[0][0]
            c = get_canonical_residue(ring)
            p_coords = (cat['coords'][0].item(), cat['coords'][1].item(), cat['coords'][2].item())
            l_coords = (ring['centroid'][0].item(), ring['centroid'][1].item(), ring['centroid'][2].item())
            t = Interaction(c, p_coords, l_coords, dist, None)

            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(I_TYPE_PI_CATION + ':' + c)


    if inters:
        return InteractionType(I_NAME_PI_CATION, list(inters.values()))
    else:
        return None


def calc_halogen_bond_interactions(protein, mol, key_inters_defs, mol_key_inters, filter_strict=False):
    protein_atoms, ligand_atoms, strict = oddt.interactions.halogenbonds(protein, mol)
    inters = {}
    for p, l, s in zip(protein_atoms, ligand_atoms, strict):
        if s or not filter_strict:
            c = get_canonical_residue(p)
            dist = spatial.distance(np.array([p['coords']]), np.array([l['coords']]))[0][0]
            p_coords = (p['coords'][0].item(), p['coords'][1].item(), p['coords'][2].item())
            l_coords = (l['coords'][0].item(), l['coords'][1].item(), l['coords'][2].item())
            t = Interaction(c, p_coords, l_coords, dist, l['id'].item())
            if c in inters:
                current_value = inters[c]
                if dist < current_value.distance:
                    inters[c] = t
            else:
                inters[c] = t
                if key_inters_defs and c in key_inters_defs:
                    mol_key_inters.append(I_TYPE_HALOGEN + ':' + c)

    if inters:
        return InteractionType(I_NAME_HALOGEN, list(inters.values()))
    else:
        return None


def get_canonical_hbond(atom):
    # print('classifying', atom['atomtype'], atom['isbackbone'], atom['isacceptor'], atom['isdonor'], atom['isdonorh'])
    res = atom['resname'] + str(atom['resnum'])
    if atom['isbackbone']:
        if atom['atomtype'] == 'N.am' or atom['atomtype'] == 'N.3':
            return res + 'BN'
        elif atom['atomtype'] == 'O.2':
            return res + 'BO'
        else:
            utils.log('Unexpected H-bond atom', res, atom['atomtype'])
    else:
        return res + 'SC'


def get_canonical_residue(atom):
    return atom['resname'] + str(atom['resnum'])


def test_read_write_json():
    with open("report.json", "r") as f:
        txt = f.read()
        data = from_json(txt)
        print('Found', len(data), 'items')
    with open("report2.json", 'w') as report:
        json.dump(data, report, cls=InteractionEncoder)
    for iset1 in data:
        for iset2 in data:
            print('Comparing', iset1.ligand_name, iset2.ligand_name)
            score, count = iset1.compare(iset2)


# def main():
#     compare_interactions(sys.argv[1], sys.argv[2])
#
#
# if __name__ == "__main__":
#     main()
