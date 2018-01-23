"""
Text Writer.
"""


def write_line(fp, line, delimit=' '):
    """
    Write given list of strings to file on a single line with an appended new line character.

    :param fp: File Object
    :param line: list or str, text to write out.
                 If type is list, joins the string elements using the specified delimiter.
    :param delimit: str, text to use to separate out line entry. Defaults to single space.
    """
    if isinstance(line, list):
        line = delimit.join(line)
    # End if
    if len(line) > 0:
        line += '\n'
        fp.write(line)
    # End if
# End write_line()


def write_multiline(fp, lines, delimit=' '):
    """
    Write given list of line strings to a file.

    :param fp: File Object
    :param lines: list, of lists with each element representing a single line
    :param delimit: str, text to use to separate line entry. Defaults to single space.
    """
    all_lines = []
    for line in lines:
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
                 'test multi string\n', 'test multi line\n', 'output that spans across two lines\n',
                 'test\tmulti\ttab-delimited\n', 'line\toutput\n']
    with open(fname, 'w') as f:
        write_line(f, 'test single string')
        write_line(f, ['test', 'list', 'single', 'line'])
        write_line(f, ['test', 'line', 'tab', 'delimited'], delimit='\t')
        write_multiline(f, [['test', 'multi', 'string']])
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
        print("{} \n should be \n {}".format(lines, pass_case))
    # End try

    os.remove(fname)
