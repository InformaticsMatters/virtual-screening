# Job: sdf-manipulate

This describes how to run the `sdf-manipulate` job from the `miscellaneous` category in the `rdkit` collection.

## What the job does

This allows to filter a SD-file based of dynamic expressions.
This uses Python dynamic expressions to filter the SD-file based on testing each record's fields using one or more
filter expressions.

## Implementation details

* Job implementation: [sdf_manip.py](/sdf_manip.py)
* Job definition: `jobs.sdf-manipulate` in [rdkit.yaml](../rdkit.yaml)

## How to run the job

### Inputs

* **Molecules**: The molecules to filter

### Options

* **Output file name**: the name of the output SD-file.
* **Filter term(s)**: one or more filter terms
* **Field types**: the data types of the files being filtered
* **Use logical OR**: combine filters using the OR logical operator.

### Outputs

The file specified by the *Output file name* option is created containing all the original data but only for the records
that pass the filter.

### Examples

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

#### 1: field names that are not valid Python variables.
The SD-file format has virtually no restrictions on the characters that can be used, but these need to be turned into
valid Python variables so that their values can be filtered. For instance a field name can contain spaces, such as
"Conformer count". To handle this all characters in the fild name that are not A-Z (upper or lower case), 0-9 or the
underscore character (_) are converted to the underscore character. So in this case, if your SD-file contains fields
named "Conformer count" then the name you must use in the --filter or --types arguments should be "Conformer_count".

#### 2: mathematical operators
A number of additional operators can also be used. These are best illustrated by examples:

"a > b "            - field a's value must be greater than field b's
"min(a, b) > 4"     - the min() function
"max(a, b) > 4"     - the max() function
"sum(a, b) > 4"     - the sum() function
"round(a, b) > 4"   - the round() function
"a < pi"            - pi comes from the Python math class whose properties and methods are available. e.g.
"floor(a) < 5.2"    - floor comes from the math class e.g. math.floor()

If additional functions are needed then please let us know, but ...

#### 3: Missing data
If a record does not contain one of the fields defined in the expression then the record fails the filter.
e.g. if you have this expression 'a < 6' and a record does not have the field 'a' then that record will fail the
filter.

#### 4: Complex expressions
Whilst most filters will probably be simple ones like 'a < 6' much more complex expressions are possible.
For instance, you could have a filter like '(a + b) > min(c, d) / 2`.

#### 5: Security
To prevent potential hacking by malicious users the filter syntax that can be used must be strictly controlled.
The functions that you can call are restricted. e.g min(a, b) is allowed, but 'exit(999)' is not.
Here we use the Python eval() function, which is a known security risk. The approach we take is twofold:
1. minimise the risk from attack. We base this mostly around this excellent article:
https://realpython.com/python-eval-function/
The code can be inspected (see [sdf_manip.py](/sdf_manip.py). If you see loopholes then please let us know.
2. In Squonk all execution takes places in a container, so, without security leaks in Docker itself,
the only harm that can be caused is to your own container. e.g. you can cause harm to yourself but not others.
You could potentially delete all your project's data, but as you already have access to that data there are many
easier ways to do so.
