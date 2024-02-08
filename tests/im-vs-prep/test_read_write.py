import rdkit_utils


def sdf_reader(file, id_col_name='_Name', recs_to_read=10):
    return rdkit_utils.SdfReader(file, id_col_name, recs_to_read=recs_to_read)


def smi_reader(file, read_header, delimiter, id_col=1):
    return rdkit_utils.SmilesReader(file, read_header, delimiter, id_col)


def smi_writer(file, delimiter):
    return rdkit_utils.SmilesWriter(file, delimiter, extra_field_names)


def test_smiles_reader_without_header():
    reader = smi_reader('data/10.smi', False, "\t")
    extra_field_names = reader.get_extra_field_names()
    assert reader.field_names is None
    assert extra_field_names is None
    assert reader.get_mol_field_name() is None
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
            assert t[3] is not None
    assert line_count == 10


def test_smiles_reader_with_header():
    reader = smi_reader('data/10+H.smi', True, "\t")
    extra_field_names = reader.get_extra_field_names()
    assert len(reader.field_names) == 2
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
            assert t[3] is not None
    assert line_count == 10


def test_read_smiles_no_header_write_smiles_no_header(tmp_path):
    reader = smi_reader('data/10.smi', False, "\t")
    out = tmp_path / 'foo.smi'
    writer = rdkit_utils.SmilesWriter(out, "\t", reader.get_extra_field_names())
    while True:
        t = reader.read()
        if not t:
            break
        mol, smi, id, props = t
        writer.write(smi, mol, id, props, [1, 2])

    writer.close()

    with open(out, 'r') as fp:
        assert len(fp.readlines()) == 10


def test_generate_header_from_smiles_with_header_no_extra_fields():
    reader = smi_reader('data/10+H.smi', True, "\t")
    headers = rdkit_utils.generate_header_values(
        reader.get_mol_field_name(), reader.field_names, 0, [])
    assert len(headers) == 2
    assert headers[0] == 'SMILES'
    assert headers[1] == 'ID'


def test_generate_header_from_smiles_with_header_2_extra_fields():
    reader = smi_reader('data/10+H.smi', True, "\t")
    headers = rdkit_utils.generate_header_values(
        reader.get_mol_field_name(), reader.field_names, 1, ['A', 'B'])
    assert len(headers) == 4
    assert headers[0] == 'SMILES'
    assert headers[1] == 'ID'
    assert headers[2] == 'A'
    assert headers[3] == 'B'


def test_generate_header_from_smiles_no_header_no_extra_fields():
    reader = smi_reader('data/10.smi', False, "\t")
    headers = rdkit_utils.generate_header_values(
        reader.get_mol_field_name(), reader.field_names, 1, [])
    assert len(headers) == 2
    assert headers[0] == 'smiles'
    assert headers[1] == 'field2'


def test_generate_header_from_smiles_no_header_2_extra_fields():
    reader = smi_reader('data/10.smi', False, "\t")
    headers = rdkit_utils.generate_header_values(
        reader.get_mol_field_name(), reader.field_names, 1, ['A', 'B'])
    assert len(headers) == 4
    assert headers[0] == 'smiles'
    assert headers[1] == 'field2'
    assert headers[2] == 'A'
    assert headers[3] == 'B'


def test_sdf_reader():
    reader = sdf_reader('data/candidates-10.sdf')
    extra_field_names = reader.get_extra_field_names()
    assert len(reader.field_names) == 5
    assert len(extra_field_names) == 5
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
            assert len(t[3]) == len(reader.field_names)
    assert molcount == 10
