import os
import subprocess
import multiprocessing
import signal, psutil
from socket import *
host = ""
port = 13000
buf = 1024
galaxy_pid = 12234
addr = (host, port)
UDPSock = socket(AF_INET, SOCK_DGRAM)
UDPSock.bind(addr)
print "Waiting to receive messages..."
while True:
    (data, addr) = UDPSock.recvfrom(buf)
    print "Received message: " + data
    if data == "exit":
        break
    elif data == "pull":
	git_pull = subprocess.Popen(["git", "pull", "origin", "fix_marshal"], stdout=subprocess.PIPE)
	output = git_pull.communicate()[0]
	print(output)
    elif data == "kill":
	for child in galaxy.children(recursive=True):
	   child.send_signal(signal.SIGINT)
    elif data == "start":
	galaxy = psutil.Popen(["/bin/sh", "/home/ABTLUS/mateus.ruivo/galaxy-dist/run.sh"], )
	#galaxy = multiprocessing.Process(target=os.system, args=('sh /home/ABTLUS/mateus.ruivo/galaxy-dist/run.sh',))
	#galaxy.start()
	print(os.getpid())
	print(galaxy.pid)
	
UDPSock.close()
os._exit(0)
