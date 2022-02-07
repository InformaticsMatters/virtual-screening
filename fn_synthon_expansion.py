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

# Cypher query that is executed:
# MATCH (fa:F2 {smiles: $smiles})
# -[:FRAG*0..%(number_hops)d]-(:F2)"  # to increase the number of hops, change '0..2' to '0..3'
# <-[e:FRAG]-(c:Mol) WHERE
# c.hac > 15 AND
# (split(e.label, '|')[1] = $synthon OR split(e.label, '|')[4] = $synthon)
# RETURN DISTINCT c"%{"number_hops":number_hops})
#

# Ruben's spec:
# The user provides the required information for execution:  a list of pairs of fragmentIds.
# Fragments could be XChem ids or could be also uploaded, as a sdf file by the user.
# List of fragment ids could be directly pasted or uploaded as file
# The user can also change the default values for querying: Number of hops, maximunSize,..
# and for filters: maxRMSD, maxEnergyGAP, overlappingFilterThr, etc.

# The system will execute the query in the Fragment network that will retrieve analogues to A containing a substructure
# of B and will do the same for analogues of B containing substructures of A.

# The system will run the filters on the molecules (smiles) retrieved by the database search.
# The last step of this step is Fragmenstein

# The system notifies run completion to the user and provides an inspect results link and download raw files link.
# Results are stored in database. Results consist of:
# - The predicted pose
# - DeltaDeltaG
# - MRMSD
# - SuCOS
# - InteractionFingerprints
# Other useful metrics could be added such as SAScore or SuCOS

# Users inspect suggested results (several result per fragmentId pair), consisting in 3D coordinates, and filters values.
# User can visualize the results and sort them by the metrics and potentially remove some of the suggested candidates
# either manually or using a thresholds

"""
Use the fragment network to find expansions of a molecule that contain a particular synthon
"""

import os, argparse, time
from rdkit import RDLogger
from neo4j import GraphDatabase
import utils
from standardize_molecule import standardize_to_noniso_smiles


RDLogger.logger().setLevel(RDLogger.ERROR)


def appendWhereClause(where_clause, prefix, value):
    if value is not None:
        if where_clause:
            where_clause += ' AND'
        where_clause += prefix + str(int(value))
    return where_clause


def process(smiles, synthon, outfile,
            neo4j_server, neo4j_username, neo4j_password,
            hac_min=None, hac_max=None,
            rac_min=None, rac_max=None,
            hops=2,
            standardize=True,
            report_hits=False
            ):

    if hops > 4:
        raise ValueError('Hops greater than 4 are not allowed')

    if standardize:
        std_smi, mol = standardize_to_noniso_smiles(smiles)
        smiles = std_smi
        std_smi, mol = standardize_to_noniso_smiles(synthon)
        synthon = std_smi
        utils.log("Standardized molecules:", smiles, synthon)

    ihops = int(hops)

    where_clause = ''
    where_clause = appendWhereClause(where_clause, ' c.hac >= ', hac_min)
    where_clause = appendWhereClause(where_clause, ' c.hac <= ', hac_max)
    where_clause = appendWhereClause(where_clause, ' c.chac >= ', rac_min)
    where_clause = appendWhereClause(where_clause, ' c.chac <= ', rac_max)

    if where_clause:
        where_clause += ' AND'

    query = "MATCH (fa:F2 {smiles: $smiles})-[:FRAG*0.." + str(ihops) + "]-(:F2)<-[e:FRAG]-(c:Mol) WHERE" +\
            where_clause +\
             " (split(e.label, '|')[1] = $synthon OR split(e.label, '|')[4] = $synthon) " +\
             "RETURN DISTINCT c"
    utils.log("QUERY:", query)

    driver = GraphDatabase.driver(neo4j_server, auth=(neo4j_username, neo4j_password))
    with driver.session() as session:
        results = session.read_transaction(run_query, query, smiles, synthon, report_hits)

    if outfile:
        with open(outfile, 'wt') as out:
            for result in results:
                out.write(result + '\n')

    return len(results)


def run_query(tx, query, smiles, synthon, report_hits):

    t0 = time.time()
    result = tx.run(query, smiles=smiles, synthon=synthon)
    t1 = time.time()
    expansions = set()
    for record in result:
        node = record['c']
        smi = node['smiles']
        if report_hits:
            utils.log(smi)
        expansions.add(smi)
    utils.log("Found", len(expansions), "molecules in", (t1 - t0), 'secs')
    return expansions


def main():
    # Example usage:
    # ./fn_synthon_expansion.py -s 'c1ccco1' -x '[Xe]C1C=CCCC1' --server "bolt://localhost:7687" \
    #    --username <username> --password <password> --hops 2

    parser = argparse.ArgumentParser(description='Fragnet synthon expansion')
    parser.add_argument('-s', '--smiles', required=True, help='Query SMILES')
    parser.add_argument('-x', '--synthon', required=True, help='Synthon to be included')
    parser.add_argument('-o', '--outfile', help='Output file')
    parser.add_argument('--report-hits', action='store_true', help='Write the results to stdout')
    parser.add_argument('--hac-min', type=int, help='The min change in heavy atom count')
    parser.add_argument('--hac-max', type=int, help='The max change in heavy atom count')
    parser.add_argument('--rac-min', type=int, help='The min change in ring atom count')
    parser.add_argument('--rac-max', type=int, help='The max change in ring atom count')
    parser.add_argument('--hops', type=int, default=2, help='The number of graph traversals (hops). Max 4')
    parser.add_argument('--no-standardize', action='store_true', help='The molecules do not need standardizing')
    parser.add_argument('--server', help='Neo4j database server')
    parser.add_argument('--username', help='Neo4j username')
    parser.add_argument('--password', help='Neo4j password')

    args = parser.parse_args()
    utils.log_dm_event("SynthonExpansion Args: ", args)

    if args.server:
        neo4j_server = args.server
    else:
        neo4j_server = os.getenv('NEO4J_SERVER')
    if not neo4j_server:
        utils.log('WARNING: no database server URL provided')

    if args.username:
        neo4j_username = args.username
    else:
        neo4j_username = os.getenv('NEO4J_USERNAME')
    if not neo4j_username:
        utils.log('WARNING: no database username provided')

    if args.password:
        neo4j_password = args.password
    else:
        neo4j_password = os.getenv('NEO4J_PASSWORD')
    if not neo4j_password:
        utils.log('WARNING: no database password provided')

    # this does the processing
    count = process(args.smiles, args.synthon, args.outfile,
                    neo4j_server, neo4j_username, neo4j_password,
                    hac_min=args.hac_min, hac_max=args.hac_max,
                    rac_min=args.rac_min, rac_max=args.rac_max,
                    hops=args.hops, standardize=not args.no_standardize,
                    report_hits=args.report_hits)

    utils.log_dm_event('Search found {} molecules.'.format(count))
    utils.log_dm_cost(count)


if __name__ == "__main__":
    main()
