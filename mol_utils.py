#!/usr/bin/env python

# Copyright 2021 Informatics Matters Ltd.
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

from rdkit import Chem

def updateChargeFlagInAtomBlock(mb):
    """
    Add data for the charges to the atom block. This data is now deprecated and should be specified using "M  CHG" lines
    but some old software such as rDock only handle the old syntax. RDKit only supports the new syntax so this method
    handles adding the old syntax as well.
    This function is based on work by Jose Manuel Gally that can be found here:
    See https://sourceforge.net/p/rdkit/mailman/message/36425493/
    """
    f="{:>10s}"*3+"{:>2}{:>4s}"+"{:>3s}"*11
    chgs = []    # list of charges
    lines = mb.split("\n")
    #if mb[0] == '' or mb[0] == "\n":
    #    del lines[0]
    CTAB = lines[3]
    atomCount = int(CTAB.split()[0])
    # parse mb line per line
    for l in lines:
        # look for M CHG property
        if l[0:6] == "M  CHG":
            records = l.split()[3:]    # M  CHG X is not needed for parsing, the info we want comes afterwards
            # record each charge into a list
            for i in range(0,len(records),2):
                idx = records[i]
                chg = records[i+1]
                chgs.append((int(idx), int(chg)))    # sort tuples by first element?
            break    # stop iterating

    # sort by idx in order to parse the molblock only once more
    chgs = sorted(chgs, key=lambda x: x[0])

    # that we have a list for the current molblock, attribute each charges
    for chg in chgs:
        i=4
        while i < 4+atomCount:    # do not read from beginning each time, rather continue parsing mb!
            # when finding the idx of the atom we want to update, extract all fields and rewrite whole sequence
            if i-3 == chg[0]:    # -4 to take into account the CTAB headers, +1 because idx begin at 1 and not 0
                fields = lines[i].split()
                x=fields[0]
                y=fields[1]
                z=fields[2]
                symb=fields[3]
                massDiff=fields[4]
                charge=fields[5]
                sp=fields[6]
                hc=fields[7]
                scb=fields[8]
                v=fields[9]
                hd=fields[10]
                nu1=fields[11]
                nu2=fields[12]
                aamn=fields[13]
                irf=fields[14]
                ecf=fields[15]
                # update charge flag
                if chg[1] == -1:
                    charge = '5'
                elif chg[1] == -2:
                    charge = '6'
                elif chg[1] == -3:
                    charge = '7'
                elif chg[1] == 1:
                    charge = '3'
                elif chg[1] == 2:
                    charge = '2'
                elif chg[1] == 3:
                    charge = '1'
                else:
                    print("ERROR! " + str(lines[0]) + "unknown charge flag: " + str(chg[1]))    # print name then go to next chg
                    break
                # update modatom block line
                lines[i] = f.format(x,y,z,symb,massDiff,charge,sp,hc,scb,v,hd,nu1,nu2,aamn,irf,ecf)
            i+=1
    upmb = "\n".join(lines)
    return(upmb)


def get_num_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    undef_cc = 0
    for cc in chiral_centers:
        if cc[1] == '?':
            undef_cc += 1
    return len(chiral_centers), undef_cc


def get_num_sp3_centres(mol):
    return sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms())

def sdf_read_mol(input):
    txt = sdf_read_block(input)
    if txt == None:
        return None
    mol = Chem.MolFromMolBlock(txt)
    print('NumProps:', mol.GetNumProps())
    return (txt, mol)


def sdf_record_gen(hnd):
    """A generator for text records fom a SD file
    """
    mol_text_tmp = ""
    while 1:
        line = hnd.readline()
        if not line:
            return
        line = line.decode("utf-8")
        mol_text_tmp += line
        if line.startswith("$$$$"):
            mol_text = mol_text_tmp
            mol_text_tmp = ""
            yield mol_text
