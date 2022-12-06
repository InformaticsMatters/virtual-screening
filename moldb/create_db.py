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
Creates the MolDB database tables.
Prior to running this the database needs to be set up. In a shell do this:

1. createdb moldb
2. psql -c 'create extension rdkit' moldb

Then set the POSTGRES_SERVER, POSTGRES_DATABASE, POSTGRES_USERNAME and POSTGRES_PASSWORD environment variables.

After creating the tables do this:

1. alter table molecule add column mol mol;
2. create index molidx on molecule using gist(mol);

TODO - clean up this creation process.
"""
from . import models

engine = models.get_engine()

models.Base.metadata.create_all(engine)

