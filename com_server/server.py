import os
import subprocess
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
	git_pull = subprocess.Popen(["git", "pull"], stdout=subprocess.PIPE)
	output = process.communicate()[0]
    elif data == "kill":
	os.kill(galaxy_pid, signal.SIGUSR1)
    elif data == "start":
    	galaxy = subprocess.Popen("sh ../../../run.sh", shell=True)
	galaxy_pid = galaxy.pid
UDPSock.close()
os._exit(0)
