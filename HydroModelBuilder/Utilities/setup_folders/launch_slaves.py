import subprocess
import sys
import os

def launch_slaves(ip_address, folders):
    """

    :param ip_address: 
    :param folders: 

    """

    cwd = os.getcwd()

    with open('slave_launcher_commands.sh', 'w') as f:
        for folder in folders:
            folder = folder.strip('\n')
            f.write('cd %s \n' %(folder))
            f.write('/home/part0075/PEST/beopest test.pst /p1 /H %s:4004 > output.txt & \n' %(ip_address))
            f.write('cd .. \n')

    #for folder in folders:
    #    folder = folder.strip('\n')    
    #    os.chdir(folder)
    #    print 'launching slave in ', os.getcwd()
    #    command = '/home/part0075/PEST/beopest test.pst /p1 /H %s:4004 > output.txt &' %(ip_address)
    #    print command   
    #    subprocess.call(command)
    #    print 'command executed'
    #    os.chdir(cwd)

if __name__ == "__main__":
    args = sys.argv
    ip_file = args[1]
    folders_file = args[2]
    with open(ip_file, 'r') as f:
        ip = f.readlines()[0]
	ip = ip.strip('\n')    

    with open(folders_file, 'r') as f:
        folders = f.readlines()

    launch_slaves(ip, folders)
    