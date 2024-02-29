import rdkit_utils


def sdf_reader(file, id_col_name='_Name', recs_to_read=10):
    return rdkit_utils.SdfReader(file, id_col_name, recs_to_read=recs_to_read)


def smi_reader(file, read_header, delimiter, id_col):
    return rdkit_utils.SmilesReader(file, read_header, delimiter, id_col, 50)


def smi_writer(file, delimiter, extra_field_names):
    return rdkit_utils.SmilesWriter(file, delimiter, extra_field_names)


def test_smiles_reader_without_header_with_id():
    reader = smi_reader('data/10.smi', False, "\t", 1)
    extra_field_names = reader.get_extra_field_names()
    assert reader.field_names == ['SMILES', 'ID']
    assert len(extra_field_names) == 0
    assert reader.get_mol_field_name() == 'SMILES'
    line_count = 0
    while True:
        t = reader.read()
        if not t:
            break
        else:
            line_count += 1
            assert t[0] is not None
            assert t[1] is not None
            assert t[2] is not None
            assert len(t[3]) == 0
    assert line_count == 10


def test_smiles_reader_without_header_no_id():
    reader = smi_reader('data/10.smi', False, "\t", None)
    extra_field_names = reader.get_extra_field_names()
    assert reader.field_names == ['SMILES', 'field2']
    assert len(extra_field_names) == 1
    assert reader.get_mol_field_name() == 'SMILES'
    line_count = 0
    while True:
        t = reader.read()
        if not t:
            break
        else:
            line_count += 1
            assert t[0] is not None
            assert t[1] is not None
            assert t[2] is None
            assert len(t[3]) == 1
    assert line_count == 10


def test_smiles_reader_with_header_with_id():
    reader = smi_reader('data/10+H.smi', True, "\t", 1)
    extra_field_names = reader.get_extra_field_names()
    assert reader.field_names == ['SMILES', 'ID']
    assert len(extra_field_names) == 0
    assert reader.get_mol_field_name() == 'SMILES'
    line_count = 0
    while True:
        t = reader.read()
        if not t:
            break
        else:
            line_count += 1
            assert t[0] is not None
            assert t[1] is not None
            assert t[2] is not None
            assert len(t[3]) == 0
    assert line_count == 10


def test_smiles_reader_with_header_no_id():
    reader = smi_reader('data/10+H.smi', True, "\t", None)
    extra_field_names = reader.get_extra_field_names()
    assert reader.field_names == ['SMILES', 'ID']
    assert len(extra_field_names) == 1
    assert reader.get_mol_field_name() == 'SMILES'
    line_count = 0
    while True:
        t = reader.read()
        if not t:
            break
        else:
            line_count += 1
            assert t[0] is not None
            assert t[1] is not None
            assert t[2] is None
            assert len(t[3]) == 1
    assert line_count == 10

# combinations
# read smi: +/- header, +/- id, +/- extra fields
#   -> write smi: +/- header, +/- omit
#   -> write sdf: +/- omit
# read sdf: +/-/_Name id +/- extra fields
#   -> write smi: +/- header, +/- omit
#   -> write sdf: +/- omit
#
# 1  -+- --
# 7  -+- +-
# 8  -+- ++
# 2  +++ +-
# 3  +++ ++
# 4  -++ +-
# 9  -++ ++
# 5  --+ +-
# 6  --+ ++
# 10 --- +-
# 11 --- ++
# 12 ++- +-
# 13 ++- ++
# 14 +-+ +-
# 15 +-+ ++


def test_read_smiles_no_header_write_smiles_no_header(tmp_path):  # 1
    # read smi: - header, + id, - extra fields    -+-
    #   -> write smi: - header, - omit            --
    reader = smi_reader('data/10.smi', False, "\t", 1)
    assert reader.field_names == ['SMILES', 'ID']
    assert reader.get_extra_field_names() == []
    out = tmp_path / 'foo.smi'
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        writer.write(smi, mol, id, props, [1, 2])

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 10
        assert len(content[0].split("\t")) == 4


def test_read_smiles_header_extra_fields_write_smiles_with_header(tmp_path):  # 2
    # read smi: + header, + id, + extra fields   +++
    #   -> write smi: + header, - omit           +-
    reader = smi_reader('data/10+H+extra.smi', True, "\t", 1)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID', 'A', 'B']
    assert reader.get_extra_field_names() == ['A', 'B']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                               reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'A', 'B', 'C', 'D']
        assert len(content[1].split("\t")) == 6


def test_read_smiles_header_extra_fields_write_smiles_with_header_omit(tmp_path):  # 3
    # read smi: + header, + id, + extra fields   +++
    #   -> write smi: + header, + omit           ++
    reader = smi_reader('data/10+H+extra.smi', True, "\t", 1)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID', 'A', 'B']
    assert reader.get_extra_field_names() == ['A', 'B']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_extra_fields_write_smiles_with_header(tmp_path):  # 4
    # read smi: - header, + id, + extra fields    -++
    #   -> write smi: + header, - omit            +-
    reader = smi_reader('data/10-extra.smi', False, "\t", 1)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID', 'field3', 'field4']
    assert reader.get_extra_field_names() == ['field3', 'field4']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'field3', 'field4', 'C', 'D']
        assert len(content[1].split("\t")) == 6


def test_read_smiles_no_id_extra_fields_write_smiles_with_header(tmp_path):  # 5
    # read smi: - header, - id, + extra fields    --+
    #   -> write smi: + header, - omit            +-
    reader = smi_reader('data/10-extra.smi', False, "\t", None)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'field2', 'field3', 'field4']
    assert reader.get_extra_field_names() == ['field2', 'field3', 'field4']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'field2', 'field3', 'field4', 'C', 'D']
        assert len(content[1].split("\t")) == 6


def test_read_smiles_no_id_extra_fields_write_smiles_with_header_omit(tmp_path):  # 6
    # read smi: - header, - id, + extra fields    --+
    #   -> write smi: + header, + omit            ++
    reader = smi_reader('data/10-extra.smi', False, "\t", None)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'field2', 'field3', 'field4']
    assert reader.get_extra_field_names() == ['field2', 'field3', 'field4']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'C', 'D']
        assert len(content[1].split("\t")) == 3


def test_read_smiles_no_header_write_smiles_with_header(tmp_path):  # 7
    # read smi: - header, + id, - extra fields    -+-
    #   -> write smi: + header, - omit            +-
    reader = smi_reader('data/10.smi', False, "\t", 1)
    assert reader.field_names == ['SMILES', 'ID']
    assert reader.get_extra_field_names() == []
    out = tmp_path / 'foo.smi'
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_no_header_write_smiles_with_header_omit(tmp_path):  # 8
    # read smi: - header, + id, - extra fields    -+-
    #   -> write smi: + header, - omit            ++
    reader = smi_reader('data/10.smi', False, "\t", 1)
    assert reader.field_names == ['SMILES', 'ID']
    assert reader.get_extra_field_names() == []
    out = tmp_path / 'foo.smi'
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_extra_fields_write_smiles_with_header_omit(tmp_path):  # 9
    # read smi: - header, + id, + extra fields    -++
    #   -> write smi: + header, + omit            ++
    reader = smi_reader('data/10-extra.smi', False, "\t", 1)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID', 'field3', 'field4']
    assert reader.get_extra_field_names() == ['field3', 'field4']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_no_header_extra_write_smiles_no_header(tmp_path):  # 10
    # read smi: - header, - id, - extra fields    ---
    #   -> write smi: - header, - omit            +-
    reader = smi_reader('data/10.smi', False, "\t", None)
    assert reader.field_names == ['SMILES', 'field2']
    assert reader.get_extra_field_names() == ['field2']
    out = tmp_path / 'foo.smi'
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'field2', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_no_header_write_smiles_no_header_omit(tmp_path):  # 11
    # read smi: - header, - id, - extra fields    ---
    #   -> write smi: - header, - omit            ++
    reader = smi_reader('data/10.smi', False, "\t", None)
    assert reader.field_names == ['SMILES', 'field2']
    assert reader.get_extra_field_names() == ['field2']
    out = tmp_path / 'foo.smi'
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'C', 'D']
        assert len(content[1].split("\t")) == 3


def test_read_smiles_header_fields_write_smiles_with_header(tmp_path):  # 12
    # read smi: + header, + id, - extra fields   ++-
    #   -> write smi: + header, - omit           +-
    reader = smi_reader('data/10+H.smi', True, "\t", 1)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID']
    assert reader.get_extra_field_names() == []
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_header_fields_write_smiles_with_header_omit(tmp_path):  # 13
    # read smi: + header, + id, - extra fields   ++-
    #   -> write smi: + header, + omit           ++
    reader = smi_reader('data/10+H.smi', True, "\t", 1)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID']
    assert reader.get_extra_field_names() == []
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is not None
        if count == 0:
            headers = rdkit_utils.generate_headers(1, 1,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'C', 'D']
        assert len(content[1].split("\t")) == 4


def test_read_smiles_header_no_id_extra_fields_write_smiles_with_header(tmp_path):  # 14
    # read smi: + header, - id, + extra fields   +-+
    #   -> write smi: + header, - omit           +-
    reader = smi_reader('data/10+H+extra.smi', True, "\t", None)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID', 'A', 'B']
    assert reader.get_extra_field_names() == ['ID', 'A', 'B']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], False)
            writer.write_header(headers)
        writer.write(smi, mol, id, props, [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'ID', 'A', 'B', 'C', 'D']
        assert len(content[1].split("\t")) == 6


def test_read_smiles_header_no_id_extra_fields_write_smiles_with_header_omit(tmp_path):  # 15
    # read smi: + header, - id, + extra fields   +-+
    #   -> write smi: + header, + omit           ++
    reader = smi_reader('data/10+H+extra.smi', True, "\t", None)
    out = tmp_path / 'foo.smi'
    assert reader.field_names == ['SMILES', 'ID', 'A', 'B']
    assert reader.get_extra_field_names() == ['ID', 'A', 'B']
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    count = 0
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        assert mol is not None
        assert smi is not None
        assert id is None
        if count == 0:
            headers = rdkit_utils.generate_headers(0, None,
                                                   reader.get_mol_field_name(), reader.field_names,
                                                   ['C', 'D'], True)
            writer.write_header(headers)
        writer.write(smi, mol, id, [], [1, 2])
        count += 1

    writer.close()

    with open(out, 'r') as fp:
        content = fp.readlines()
        assert len(content) == 11
        assert content[0].strip().split("\t") == ['SMILES', 'C', 'D']
        assert len(content[1].split("\t")) == 3


def test_generate_header_from_smiles_with_header_no_extra_fields():
    reader = smi_reader('data/10+H.smi', True, "\t", 1)
    headers = rdkit_utils.generate_headers(1, 1,
        reader.get_mol_field_name(), reader.field_names, [], False)
    assert len(headers) == 2
    assert headers[0] == 'SMILES'
    assert headers[1] == 'ID'


def test_generate_header_from_smiles_with_header_2_extra_fields():
    reader = smi_reader('data/10+H.smi', True, "\t", 1)
    headers = rdkit_utils.generate_headers(1, 1,
        reader.get_mol_field_name(), reader.field_names, ['A', 'B'], False)
    assert len(headers) == 4
    assert headers[0] == 'SMILES'
    assert headers[1] == 'ID'
    assert headers[2] == 'A'
    assert headers[3] == 'B'


def test_generate_header_from_smiles_no_header_no_extra_fields():
    reader = smi_reader('data/10.smi', False, "\t", 1)
    headers = rdkit_utils.generate_headers(1, 1,
        reader.get_mol_field_name(), reader.field_names, [], False)
    assert len(headers) == 2
    assert headers[0] == 'SMILES'
    assert headers[1] == 'ID'


def test_generate_header_from_smiles_no_header_2_extra_fields():
    reader = smi_reader('data/10.smi', False, "\t", 1)
    headers = rdkit_utils.generate_headers(1, 1,
        reader.get_mol_field_name(), reader.field_names, ['A', 'B'], False)
    assert len(headers) == 4
    assert headers[0] == 'SMILES'
    assert headers[1] == 'ID'
    assert headers[2] == 'A'
    assert headers[3] == 'B'


def test_sdf_reader():
    reader = sdf_reader('data/candidates-10.sdf')
    extra_field_names = reader.get_extra_field_names()
    assert len(reader.field_names) == 6
    assert len(extra_field_names) == 6
    assert reader.get_mol_field_name() is None
    molcount = 0
    while True:
        t = reader.read()
        if not t:
            break
        else:
            molcount += 1
            assert t[0] is not None
            assert t[1] is not None
            assert t[2] is not None
            assert t[3] is not None
            assert len(t[3]) == 7
    assert molcount == 10
