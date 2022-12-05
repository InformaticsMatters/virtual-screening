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

import os

from sqlalchemy import Column, PrimaryKeyConstraint, ForeignKey, Index, UniqueConstraint, func
from sqlalchemy import MetaData, Table, Column, Integer, BigInteger, SmallInteger, String, DateTime, Text, Float, CHAR
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, declarative_mixin, relationship, sessionmaker

import sqlalchemy

import utils

utils.log('sqlalchemy version:', sqlalchemy.__version__)


Base = declarative_base()

@declarative_mixin
class TimestampMixin:
    created_at = Column(DateTime, default=func.now())


class Molecule(TimestampMixin, Base):

    __tablename__ = 'molecule'

    id = Column(BigInteger, primary_key=True)

    smiles = Column(Text(), nullable=False, unique=True)

    hac = Column(SmallInteger(), index=True)
    rot_bonds = Column(SmallInteger())
    ring_count = Column(SmallInteger())
    aromatic_ring_count = Column(SmallInteger())
    chiral_centres = Column(SmallInteger())
    undef_chiral_centres = Column(SmallInteger())
    num_sp3 = Column(SmallInteger())
    logp = Column(Float())
    tpsa = Column(Float())

    def __repr__(self):
        return "<Molecule(id=%s, smiles='%s', created_at='%s')>" \
               % (self.id, self.smiles, self.created_at)


class Enumeration(TimestampMixin, Base):

    __tablename__ = 'enumeration'

    id = Column(BigInteger, primary_key=True)

    molecule_id = Column(BigInteger, ForeignKey('molecule.id', ondelete="CASCADE"), nullable=False, index=True)
    code = Column(CHAR(1), nullable=False)
    smiles = Column(Text, nullable=False)
    coords = Column(Text, nullable=False)

    def __repr__(self):
        return "<Enumeration(id=%s, molecule_id='%s', code='%s', smiles='%s', created_at='%s')>"\
               % (self.id, self.molecule_id, self.code, self.smiles, self.created_at)


class Conformer(TimestampMixin, Base):

    __tablename__ = 'conformer'

    id = Column(BigInteger, primary_key=True)

    enumeration_id = Column(BigInteger, ForeignKey('enumeration.id', ondelete="CASCADE"), nullable=False, index=True)
    coords = Column(Text, nullable=False)
    energy = Column(Float)
    energy_delta = Column(Float)

    def __repr__(self):
        return "<Conformer(id=%s, enumeration_id='%s', coords='%s', created_at='%s')>" \
               % (self.id, self.enumeration_id, self.coords, self.created_at)


class Supply(TimestampMixin, Base):

    __tablename__ = 'supply'

    id = Column(BigInteger, primary_key=True)

    smiles = Column(Text(), nullable=False)
    code = Column(Text(), nullable=False)
    molecule_id = Column(BigInteger, ForeignKey('molecule.id', ondelete="CASCADE"), nullable=False, index=True)
    file_id = Column(Integer, ForeignKey('file.id', ondelete="CASCADE"), nullable=False, index=True)

    file = relationship("File", back_populates="supplies")

    def __repr__(self):
        return "<Supply(id=%s, smiles='%s', code='%s', created_at='%s')>" % (self.id, self.smiles, self.code, self.created_at)



class File(TimestampMixin, Base):

    __tablename__ = 'file'

    id = Column(Integer, primary_key=True)
    name = Column(Text(), nullable=False)
    library_id = Column(Integer, ForeignKey('library.id', ondelete="CASCADE"), nullable=False, index=True)
    load_table = Column(Text(), nullable=False)

    library = relationship("Library", back_populates="files")
    supplies = relationship("Supply", back_populates="file", cascade="all, delete-orphan")

    __table_args__ = (
        UniqueConstraint(name, library_id, name='uq_file'),
    )

    def __repr__(self):
        return "<File(id=%s, name='%s', created_at='%s', library_id='%s', load_table='%s')>" % (self.id, self.name, self.created_at, self.library_id, self.load_table)


class Library(TimestampMixin, Base):

    __tablename__ = 'library'

    id = Column(Integer, primary_key=True)
    name = Column(Text(), nullable=False, unique=True)

    files = relationship("File", back_populates="library", cascade="all, delete-orphan")

    def __repr__(self):
        return "<Library(id=%s, name='%s', created_at='%s')>" % (self.id, self.name, self.created_at)


_engine = None

def get_engine(echo=True):

    utils.log('Using database ' + gen_url(obscure_password=True))

    global _engine
    if _engine is None:
        url = gen_url(obscure_password=False)
        _engine = create_engine(url, echo=echo, future=True)
    return _engine


def gen_url(obscure_password=True):
    pg_server = os.getenv('POSTGRES_SERVER', default='localhost')
    pg_database = os.getenv('POSTGRES_DATABASE', default='postgres')
    pg_username = os.getenv('POSTGRES_USERNAME', default='postgres')
    pg_password = os.getenv('POSTGRES_PASSWORD', default='squonk')

    templ = 'postgresql://{}:{}@{}/{}'
    if obscure_password:
        url = templ.format(pg_username, '********', pg_server, pg_database)
    else:
        url = templ.format(pg_username, pg_password, pg_server, pg_database)
    return url