"""Utilities to read and write from/to an SDF file

    Note that this was copied from the sdf-format-support repo. It was copied
    here for speed but this and convertFile SHOULD be in a file accessible
    by both virtual-screening and the formatter (and possibly others).

An SDF file consists of a series of records separated by a "$$$$',
where each record consists of:
1. A molecule block (from the start of the record up to and including the 'END' line. The first
line of this block may be a header line that s definable by the user. Example:

O=C(CSCc1ccc(Cl)s1)N1CCC(O)CC1
     RDKit          3D

 18 19  0  0  0  0  0  0  0  0999 V2000
    8.7102   -1.3539   24.2760 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.4334   -2.1203   23.6716 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3260   -1.7920   22.4941 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.5607   -0.5667   21.3699 S   0  0  0  0  0  0  0  0  0  0  0  0
    7.9641   -1.3976   21.0216 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1007   -0.5241   20.1671 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7930   -0.1276   20.3932 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2841    0.6934   19.3422 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2234    0.8796   18.3624 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    1.8209   16.9402 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    7.6812    0.0795   18.6678 S   0  0  0  0  0  0  0  0  0  0  0  0
    9.5928   -3.4405   24.2306 N   0  0  0  0  0  0  0  0  0  0  0  0
   10.8197   -3.4856   25.0609 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0016   -4.9279   25.4571 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.9315   -5.2800   26.4615 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.3887   -4.7677   27.7090 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5793   -4.6419   26.1747 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3826   -4.0949   24.7695 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  2 12  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
  9 11  1  0
 11  6  1  0
 12 13  1  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 15 17  1  0
 17 18  1  0
 18 12  1  0
M  END

2. A series of formatted parameter lines. Example:

>  <TransFSScore>  (1)
0.115601

where (1) is a record number.

"""

import re
import uuid
from typing import Dict, List

_prop_pattern: str = '^>\\s+<(.*)>'
_mol_separator: str = 'M  END'
_rec_separator: str = '$$$$'


# Supporting functions for Json
def is_valid_uuid(value: str):
    """"
    Checks whether ths given value is a UUID
    """

    try:
        uuid.UUID(str(value))
        return True
    except ValueError:
        return False


def sdf_add_property(properties: Dict[str, str], propname: str, propvalue: List[str]):
    """"
    Adds a property (name, value) pair to the properties dictionary
    """

    if propvalue[len(propvalue) - 1]:
        # property block should end with an empty line, but some SDFs are buggy
        properties[propname] = '\n'.join(propvalue)
    else:
        properties[propname] = '\n'.join(propvalue[:-1])

    return properties


# Supporting function for file
def sdf_get_next_record(sdf_file):
    """Gets the next sdf record from the supplied file.

    :param sdf_file: The handle of the sdf file object to process
    :returns:
    """
    molecule_block = []
    molecule_name = ''
    properties: Dict[str, str] = {}

    # Loop through file to get molecule block
    while True:

        # Get next line from file
        line = sdf_file.readline()
        # if line is empty end of file is reached
        if not line:
            return molecule_block, molecule_name, properties

        # remove newline chars
        text = line.rstrip('\n\r')
        molecule_block.append(text)
        # 'M  END' signifies the end of the molblock. The properties will follow
        if text == _mol_separator:
            break

    # molecule name is first line of molecule block if present
    if molecule_block[0]:
        molecule_name = molecule_block[0]

    # Add a line feed.
    molecule_block = '\n'.join(molecule_block)

    # Loop through properties
    property_name = None
    propvalue = []

    while True:
        line = sdf_file.readline()
        text = line.strip()

        # '$$$$' signifies the end of the record
        if text == _rec_separator:
            if property_name:
                properties = sdf_add_property(properties, property_name, propvalue)
            return molecule_block, molecule_name, properties

        # Check if text is a formatted sdf property name
        property_found = re.match(_prop_pattern, text)
        if property_found:
            # For each sdf property the value will follow on the next line
            if property_name:
                properties = sdf_add_property(properties, property_name, propvalue)
            property_name = property_found.group(1)
            propvalue = []
        else:
            # For each sdf property the value will follow on the next line
            propvalue.append(text)


def sdf_write_record(output_sdf_file, molecule_block, molecule_name, properties, rec_nr):
    """Write a formatted SDF record to the supplied SDF file handle.

    :param output_sdf_file: The handle of the sdf file object to process
    :param molecule_block: Formatter sdf molecule block
    :param molecule_name: Molecule name replaces the first line of the molecule block if
         present.
    :param properties: List of parameters to write to the file.
    :param rec_nr: Record number to write to SDF file
    :returns:
    """
    if molecule_name:
        # Replace the characters to the first newline with the uuid
        rest_molecule_block = molecule_block.split('\n', 1)[1]
        molecule_block = molecule_name+'\n'+rest_molecule_block

    output_sdf_file.write(molecule_block + '\n')

    prop_line = '>  <{name}>  ({rec_nr})' + '\n'
    for name, value in properties.items():
        output_sdf_file.write(prop_line.format(name=name, rec_nr=rec_nr))
        output_sdf_file.write(value + '\n\n')

    # record separator should be on a separate line
    output_sdf_file.write(_rec_separator+'\n')
