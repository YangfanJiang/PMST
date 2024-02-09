import subprocess
import configparser

config = configparser.ConfigParser()
param_f = "parameters.ini"
config.read(param_f)
process_num = config['Process_Num']['process_num']

cmd = 'python -m scoop -n '+process_num+' geo11.py'

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
for l in out.decode("utf-8").split('\r\n'):
    print(l)