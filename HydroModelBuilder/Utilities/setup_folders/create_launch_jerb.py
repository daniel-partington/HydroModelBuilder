import subprocess
import sys

nprocessors = 12

with open('jerb', 'w') as f:
        f.write('#$ -cwd \n')
        f.write('#$ -pe threaded %i \n' %(nprocessors))
        f.write('#$ -m be -M part0075@flinders.edu.au \n')
        f.write('#$ -q 1g.q \n')
        f.write('\n')
        f.write('/home/part0075/PEST/beopest test.pst /p1 /H :4004 > output.txt & \n')
        f.write("hostname -I | cut -f1 -d' ' > ip.txt \n")
	f.write('\n')
	f.write('python launch_slaves.py ip.txt folders.txt \n')
	f.write('chmod 777 slave_launcher_commands.sh \n')
	f.write('./slave_launcher_commands.sh \n')
        f.write('\n')
        f.write('wait \n')
        f.write('exit \n')

subprocess.call("qsub jerb" , shell=True)

# if __name__ == "__main__":
    