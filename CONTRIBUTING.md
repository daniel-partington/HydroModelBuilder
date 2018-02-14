# Contribution Guide

## Code Style

Follow PEP8 guidelines with the help of a linter, [preferably 'autopep8'](https://pypi.python.org/pypi/autopep8),
with the following exceptions:

* Max line length set to 120 characters

Plugins for common IDEs/editors.

* [https://github.com/spyder-ide/spyder-autopep8]
* [https://atom.io/packages/python-autopep8]

Install autopep8 via pip or conda

`pip install autopep8`

## Docstrings

Docstrings should follow the sphinxdoc style and are required where the function has one or more of the following:

* has input parameters
* more than 3 lines long
* returns a variable of a type not immediately determinable
* raises an exception

Example docstring:

```
"""This is a reST style docstring, modified to include type information

Uses reST style markup for sphinx to generate documentation.
For example: an equation is written like :math:`e = mc^2`

Test snippets/examples can also be included using triple arrowheads

>>> p1 = "test"
>>> p2 = ["test", "list"]
>>> this_function(p1, p2)
["the", "expected", "output", "from", "function", "call"]

:param param1: str, this is the first param
:param param2: list[str], give type of elements if they are uniform. If they are not, `object` is okay but describe what the expected input is.
:returns: this is a description of what is returned
:raises KeyError: raises an exception when...
"""
```
