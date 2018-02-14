# Contribution Guide

## Code Style

Follow PEP8 guidelines with the help of a linter, preferably ['pycodestyle'](http://pycodestyle.readthedocs.io/en/latest/),
with the following exceptions:

* Max line length set to 120 characters
* Use spaces for indentation (4 spaces) instead of tabs as is preferred.

Linter Plugins for commonly used IDEs/editors.

* https://github.com/spyder-ide/spyder-autopep8
* https://github.com/AtomLinter/linter-pycodestyle

Linters can be installed via pip or conda

`pip install autopep8`

`pip install pycodestyle`

## Docstrings

Docstrings should follow the sphinx style (which uses [reStructuredText](http://www.sphinx-doc.org/en/master/rest.html))
and are required where the function has one or more of the following:

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
:param param2: list[str], give type of elements if they are uniform.
               If they are not, `object` is okay but describe what the expected input is.

:returns: list[str], this is a description of what is returned
:raises KeyError: raises an exception when...
"""
```

Docstrings can be a single line if the function does not have any parameters, returns, and raises no exceptions:

```
"""A single line docstring"""
```
