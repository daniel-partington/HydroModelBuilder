"""
Text Writer.
"""


def write_line(fp, line, delimit=' '):
    """Write given list of strings to file on a single line with an appended new line character.

    :param fp: File Object
    :param line: list or str, text to write out.
                 If type is list, converts all elements to string and joins them using the specified delimiter.
    :param delimit: str, text to use to separate out line entry. (Default value = ' ')

    """

    if isinstance(line, list):
        line = [str(el) for el in line]
        line = delimit.join(line)
    # End if

    if len(line) > 0:
        line += '\n'
        fp.write(line)
    # End if
# End write_line()


def write_multiline(fp, lines, delimit=' '):
    """Write given list of line strings to a file.

    :param fp: File Object
    :param lines: list[list], with each list element representing a single line
    :param delimit: str, text to use to separate line entry. (Default value = ' ')
    """

    all_lines = []
    for line in lines:
        line = [str(el) for el in line]
        str_to_write = delimit.join(line) + '\n'
        all_lines.append(str_to_write)
    # End for
    if len(all_lines) > 0:
        fp.writelines(all_lines)
# End write_multiline()


if __name__ == '__main__':
    import os
    fname = 'test_textwriter.txt'
    pass_case = ['test single string\n', 'test list single line\n', 'test\tline\ttab\tdelimited\n',
                 '1000 1000 50 50.0 1000000.0\n',
                 'test multi string\n',
                 '1000 9999.9 100000.0\n',
                 'test multi line\n', 'output that spans across two lines\n',
                 'test\tmulti\ttab-delimited\n', 'line\toutput\n']

    with open(fname, 'w') as f:
        write_line(f, 'test single string')
        write_line(f, ['test', 'list', 'single', 'line'])
        write_line(f, ['test', 'line', 'tab', 'delimited'], delimit='\t')
        write_line(f, ['1000', 1000, 50, 50.0, 1e6])

        write_multiline(f, [['test', 'multi', 'string']])
        write_multiline(f, [[1000, 9999.9, 1e5]])
        write_multiline(f, [['test', 'multi', 'line'], ['output', 'that', 'spans', 'across', 'two lines']])
        write_multiline(f, [['test', 'multi', 'tab-delimited'], ['line', 'output']], delimit='\t')

    try:
        with open(fname, 'r') as f:
            lines = f.readlines()
            assert lines == pass_case, "Test Failed!"
        # End with

        print("Test passed")
    except AssertionError as e:
        print(e)
        print(("{} \n should be \n {}".format(lines, pass_case)))
    # End try

    os.remove(fname)
