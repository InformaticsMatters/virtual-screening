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


# IMPORTANT: this currently has to run using Python2 as that's all the rdock container has.

import argparse, os, subprocess, time, datetime
from jinja2 import Template

template = Template('''RBT_PARAMETER_FILE_V1.00
TITLE rDockdocking

RECEPTOR_FILE {{ receptor }}
RECEPTOR_FLEX 3.0

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {{ ligand }}
    RADIUS 3.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION

''')

# these 2 methods are duplicated from utils.py because of the pytnon2 issue
def log(msg):
    msg_time = datetime.datetime.utcnow().replace(microsecond=0)
    print('%s # INFO -EVENT- %s' % (msg_time, msg))

def expand_path(path):
    """
    Create any necessary directories to ensure that the file path is valid
    
    :param path: a filename or directory that might or not exist
    """
    head_tail = os.path.split(path)
    if head_tail[0]:
        if not os.path.isdir(head_tail[0]):
            log('Creating directories for', head_tail[0])
            os.makedirs(head_tail[0], exist_ok=True)
            

def execute(receptor, ligand, output):

    prmfilename = output + '.prm'
    #pharmafile = output + '.restr'
    content = template.render(receptor=receptor, ligand=ligand)
    
    # generate an empty pharmacophore restraint file 
    #with open(pharmafile, 'w') as pharma:
    #    pass
    
    # write the rdock configuration file (the .prm file)
    log("Generating prm file")
    expand_path(prmfilename)
    with open(prmfilename, 'w') as prmfile:
        prmfile.write(content)
        
    # run the rbcavity program to create the active site definition files
    # they will have base names of the output variable     
    cmd1 = ['rbcavity', '-was', '-d', '-r', prmfilename]
    log('Generating cavity definition: ' + ' '.join(cmd1))
    exit_code = subprocess.call(cmd1)
    
    if exit_code:
        log("Failed to generate cavity. Exit code: " + str(exit_code))
        raise Exception("Failed to generate cavity")
    log("rbcavity completed successfully")


def main():

    # Example:
    #   python3 prepare_rdock.py --receptor receptor.mol2 --ligand ligand.mol --output docking

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare rDock docking')
    parser.add_argument('-r', '--receptor', required=True, help='Prepared receptor in MOL2 format')
    parser.add_argument('-l', '--ligand', required=True, help='Ligand in Molfile format')
    parser.add_argument('-o', '--output', required=True, default='docking', help='Base name for outputs')
   

    args = parser.parse_args()
    log("prepare_rdock.py: ", args)
    
    t0 = time.time()
    execute(args.receptor, args.ligand, args.output)
    t1 = time.time()
    
    log('Generated files ' + args.output + '.prm and ' + args.output + '.as' + ' in ' + str(t1 - t0) + ' seconds') 
    
    
if __name__ == "__main__":
    main()
