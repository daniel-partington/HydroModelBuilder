" Utility to convert Unix to DOS format "

import sys

def unix2dos(path):
    """

    :param path: 

    """
    text = open(path, "U").read() 
    text = text.replace("\n", "\r\n") 
    open(path, "wb").write(text) 

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except:
        print("unix2dos requires filename as argument")

    unix2dos(path)        