
import os
import shutil

def createFolders(num, inlocation=None, prefix=None):
    """
    Function to create multiple folders with unique id's which by default 
    are just integers, but can also have a specified prefix. All folders are 
    created in the location 'inlocation' if specified otherwise it is just
    in the working directory
    
    :param name: num: number of folder to create which will be indexed 0, 1, 2 ...
    :param name: inlocation: string specifying folder location to create new folders
    :param name: prefix: optional string to place at start of folder names
    """

    if inlocation == None:
        inlocation = os.getcwd()

    if prefix == None:
        prefix = ''

    print 'Creating %i folders' %(num)   

    folder_list = []
    
    # Write folder list to a text file
    with open(inlocation + os.path.sep + 'folders.txt', 'w') as f:
        for i in range(num):
            folder_name = inlocation + os.path.sep + prefix + str(i)             
            try:        
                os.mkdir(folder_name)
                folder_list += folder_name 
                f.write(folder_name+'\n')
            except:
                print 'Could not create: ', folder_name

        
def cleanupFolders(num, inlocation=None, prefix=None):
    """
    Function to delete multiple folders created with createFolders() that contain
    unique id's as specified in createFolders()
    
    :param name: num: number of folder to delete which are indexed 0, 1, 2 ...
    :param name: inlocation: string specifying folder location to delete folders
    :param name: prefix: optional string specifying prefix string
    """

    if inlocation == None:
        inlocation = os.getcwd()

    if prefix == None:
        prefix = ''

    print 'Deleting %i folders' %(num)   
    
    for i in range(num):
        folder_name = inlocation + os.path.sep + prefix + str(i)             
        try:
            shutil.rmtree(folder_name)
        except:
            print 'Unable to remove: ', folder_name


def cleanupFoldersByList(folder_list):
    """
    Function to delete folders from a list of folders of from a text file 
    containing full path to folder strings, with one on each line
    
    :param name: folder_list: a list of folders or string containing filename with folders
    """

    if type(folder_list) == list:
        folder_list = folder_list
    elif type(folder_list) == str:
        filename = folder_list        
        with open(filename, 'r') as f:
            folder_list = f.readlines()
        
    print 'Deleting %i folders' %(len(folder_list))
    
    for folder_name in folder_list:
        try:
            shutil.rmtree(folder_name.strip('\n'))
        except:
            print 'Unable to remove: ', folder_name

# From http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
def copytree(src, dst, symlinks=False, ignore=None):
    """
    Function to copy contents of one directory to another
    :param name: src: source directory
    :param name: dst: destination directory
    :param name: symlinks: ??
    :param name: ignore: ??
    
    """
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def setupModelDirectories(num, source_directory, inlocation=None, prefix=None ):
    """ 
    Function to setup 'num' folders to run models based on files in 'source_directory'
    in location 'inlocation' and allowing for additional prefix in folder names

    :param name: num: number of folder to create which will be indexed 0, 1, 2 ...
    :param name: source_directory: string indicating folder where model resides
    :param name: inlocation: string specifying folder location to create new folders
    :param name: prefix: optional string to place at start of folder names
    """    
    # create the folders
    createFolders(num, inlocation=inlocation, prefix=prefix)
    if inlocation == None:
        inlocation = '.'

    if prefix == None:
        prefix = ''

    for i in range(num):
        folder = inlocation + os.path.sep + prefix + str(i)
        copytree(source_directory, folder)
    

if __name__ == '__main__':
    createFolders(3)
    cleanupFoldersByList('folders.txt')
