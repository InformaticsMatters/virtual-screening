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
Filter SD-file using Python dynamic expressions.

Example:
   python -m sdf_manip -i data/mols.sdf -o foo.sdf -f 'Energy < 1' 'o3da_score < 80' \
     --types "Energy:float, o3da_score:float"

Parameters:
    -i --input The input SD-file
    -o --output The output SD-file
    -f --filter The filter(s) to apply
    -t --types The field types (int, float, str)

This uses Python dynamic expressions to filter the SD-file based on testing each record's fields using one or more
filter expressions.

Simple example:
Your SD-file contains records that have a field named 'quantity' that contains decimal numbers.
You want to filter those records so that only ones with a quantity greater than 9.3 are included in the output.
The filter terms could look like this: --filter 'quantity > 9.3' and the type term like this: '--types quantity:float'.
What happens here is that the SD-file fields are exposed as 'variables' of the same name (quantity) and those values
are then filtered using the expression 'quantity > 9.3'. If the that expression evaluates to True then the record is
written to the output.
Types that are supported are int (integer), float (32 bit floating point number) and str (string or text). These are
Python data types. Specifying str is not really necessary as that is the default if you don't specify anything.

Multiple filter terms can be specified. The filter argument can take multiple arguments and the types argument can
specify the types of multiple fields. For instance consider these arguments:
    --filter "quantity < 9.3" "count > 2.4" '--types "quantity:float, count:int"'
Note that filtering uses two fields, quantity being treated as a floating point number and count as an integer.
Note the syntax: the --filter argument takes 1 or more filter expressions whiles the --types argument has a single value
whose terms are separated by a comma.

In many cases multiple filters can be combined into a single filter. e.g. the two terms "quantity < 9.3" and
"count > 2.4" can be combined into one term as "quantity < 9.3 and count > 2.4". Note the "and" keyword used to
combine the terms into a single expression.

And on the subject of how to logically combine expressions, if you do have multiple filter terms then by default they
are combined using the AND logical operator, but you can change this to OR using the --logic-or argument.

That example, though accurate, hides a few aspects.

1: field names that are not valid Python variables.
The SD-file format has virtually no restrictions on the characters that can be used, but these need to be turned into
valid Python variables so that their values can be filtered. For instance a field name can contain spaces, such as
"Conformer count". To handle this all characters in the fild name that are not A-Z (upper or lower case), 0-9 or the
underscore character (_) are converted to the underscore character. So in this case, if your SD-file contains fields
named "Conformer count" then the name you must use in the --filter or --types arguments should be "Conformer_count".

2. mathematical operators
A number of additional operators can also be used. These are best illustrated by examples:

"a > b "            - field a's value must be greater than field b's
"min(a, b) > 4"     - the min() function
"max(a, b) > 4"     - the max() function
"sum(a, b) > 4"     - the sum() function
"round(a, b) > 4"   - the round() function
"a < pi"            - pi comes from the Python math class whose properties and methods are available. e.g.
"floor(a) < 5.2"    - floor comes from the math class e.g. math.floor()

If additional functions are needed then please let us know, but ...

3. Missing data
If a record does not contain one of the fields defined in the expression then the record fails the filter.
e.g. if you have this expression 'a < 6' and a record does not have the field 'a' then that record will fail the
filter.

4. Complex expressions
Whilst most filters will probably be simple ones like 'a < 6' much more complex expressions are possible.
For instance, you could have a filter like '(a + b) > min(c, d) / 2`.

5. Security
To prevent potential hacking by malicious users the filter syntax that can be used must be strictly controlled.
The functions that you can call are restricted. e.g min(a, b) is allowed, but 'exit(999)' is not.
Here we use the Python eval() function, which is a known security risk. The approach we take is twofold:
    1. minimise the risk from attack. We base this mostly around this excellent article:
    https://realpython.com/python-eval-function/
    The code below can be inspected. If you see loopholes then please let us know.
    2. In Squonk all execution takes places in a container, so, without security leaks in Docker itself,
    the only harm that can be caused is to your own container. e.g. you can cause harm to yourself but not others.
    You could potentially delete all your project's data, but as you already have access to that data there are many
    easier ways to do so.

"""

import re, math, argparse, sys

from rdkit import Chem

import utils
from dm_job_utilities.dm_log import DmLog


ALLOWED_NAMES = {'min': min, 'max': max, 'sum': sum, 'round': round}
for k, v in math.__dict__.items():
    ALLOWED_NAMES[k] = v

ALLOWED_TYPES = {'str': str, 'int': int, 'float': float}


class Record(object):
    def __init__(self, adict, dtypes):
        d = {}
        for k, v in adict.items():
            repl = re.sub(r"\W", "_", k)
            if repl in dtypes:
                v = dtypes[repl](v)
            d[repl] = v

        self.values = d

    def test(self, term):
        return eval(term, {"__builtins__": ALLOWED_NAMES}, self.values)


def execute(input, output, filters, types, logical_or=False):

    DmLog.emit_event('Filters:', filters)
    DmLog.emit_event('Types:', types)

    num_pass = 0
    num_fail = 0
    num_errors = 0
    with open(output, 'wt') as writer:
        with Chem.SDMolSupplier(input) as supplr:
            count = -1
            for mol in supplr:
                count += 1
                d = mol.GetPropsAsDict()
                try:
                    r = Record(d, types)
                except:
                    utils.log('Failed to handle record', count)
                    num_errors += 1
                    continue
                b = []
                for filt in filters:
                    result = False
                    try:
                        result = r.test(filt)
                    except NameError as e:
                        # Also thrown when using a missing function name
                        ex_type, ex_value, ex_traceback = sys.exc_info()
                        utils.log('Error in record {}: {}'.format(count, ex_value))
                        num_errors += 1
                        b.append(False)
                    b.append(result)
                if logical_or:
                    fail = True not in b
                else:
                    fail = False in b

                if fail:
                    num_fail += 1
                else:
                    num_pass += 1
                    writer.write(supplr.GetItemText(count))

    DmLog.emit_event('Stats:', num_pass, 'pass,', num_fail, 'fail,', num_errors, 'errors')
    DmLog.emit_cost(count + 1)


def _gen_types_map(types_str):
    tokens = types_str.split(',')
    d = {}
    for token in tokens:
        parts = token.split(':')
        if not len(parts) == 2:
            raise ValueError('Unexpected type specification: ' + token)
        name = parts[0].strip()
        type_str = parts[1].strip()
        if type_str not in ALLOWED_TYPES:
            raise ValueError('Unexpected type: ' + type_str)
        d[name] = ALLOWED_TYPES[type_str]
    return d


def main():

    # Example:
    #   python -m sdf_manip -i data/mols.sdf -o foo.sdf -f 'Energy < 1' 'o3da_score < 80' \
    #     --types "Energy:float, o3da_score:float"

    # ----- command line args definitions ---------------------------------------------------------

    parser = argparse.ArgumentParser(description='SDF manipulate')
    parser.add_argument('-i', '--input', required=True,  help="Input file (.sdf)")
    parser.add_argument('-o', '--output', required=True,  help="Output file (.sdf)")
    parser.add_argument('-f', '--filter', nargs='+', help="Filter specification(s)")
    parser.add_argument('-t', '--types', help="Variable types")
    parser.add_argument('--logic-or', action='store_true', help="Combine filters using OR rather than AND")

    args = parser.parse_args()
    DmLog.emit_event("sdf manip: ", args)

    if args.types:
        types_dict = _gen_types_map(args.types)
    else:
        types_dict = {}

    execute(args.input, args.output, args.filter, types_dict, logical_or=args.logic_or)


if __name__ == "__main__":
    main()
