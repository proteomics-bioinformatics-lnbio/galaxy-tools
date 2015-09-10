import os, time
from socket import *
#host = "10.0.2.185" # set to IP address of target computer
host = "127.0.0.1"
port = 13000
addr = (host, port)
UDPSock = socket(AF_INET, SOCK_DGRAM)
while True:
    data = raw_input("Enter message to send or type 'exit': ")
    UDPSock.sendto(data, addr)
    time.sleep(5.0)
    if data == "exit":
        break
UDPSock.close()
os._exit(0)
