
import subprocess

num = 12

with open('jerb', 'w') as f:
        f.write('#$ -cwd \n')
        f.write('#$ -pe threaded %i \n' %(num))
        f.write('#$ -m be -M part0075@flinders.edu.au \n')
        f.write('#$ -q 1g.q \n')
        f.write('\n')
        f.write('../beopest test.pst /p1 /H :4004 > output.txt & \n')
        f.write('hostname -I | cut -f1 -d' ' > ip.txt')
	f.write('\n')
	f.write('python launch_slaves.py ip.txt slave_folders.txt')
        f.write('\n')
        f.write('wait \n')
        f.write('exit \n')

subprocess.call("qsub jerb" , shell=True)