#!/usr/bin/python

from scipy import fromfile
from pylab import *
from sys import argv
g=loadtxt('zgrid.out')
#a=fromfile('op.rho.0000005000.bin',float32)
#a=fromfile('op.vor.0000005000.bin',float32)
a=fromfile(argv[1],float32)
print a.shape
NI,NJ,NK = 48,32,16
NI,NJ,NK =48,48,24 
NI,NJ,NK =96,96,32
NI,NJ,NK =96,96,24
NI,NJ,NK =96,96,64
NI,NJ,NK =48,96,32
NI,NJ,NK = 96,192,32
NI,NJ,NK =48,48,32
NI,NJ,NK = 96,192,64
NI,NJ,NK = 96,192,32
z1 = g[0:NK+2,1]
print z1.shape
print g.shape
y,z = meshgrid(arange(NJ+2),z1)
print y.shape

b=a.reshape(NK+2,NJ+2,NI+2)

print 'b.shape',b.shape,'b.sum',b.sum()
figure(1)
contour(b[-8,:,:])
colorbar()
title('xy')
xlabel('x')
ylabel('y')
savefig(argv[1]+'xy.png')
figure(2)
contour(y,z,b[:,:,1])
colorbar()
xlabel('y')
ylabel('z')
title('yz')
savefig(argv[1]+'zy.png')
figure(3)
subplot(2,2,1)
plot(b[:,NJ/2,NI/2].flatten(),z[:,0])
subplot(2,2,2)
plot(b[-1,:,NI/2].flatten())
subplot(2,2,3)
plot(-10./1025.*diff(b[:,NJ/2,NI/2].flatten())/diff(z1),z1[1::])
xlabel('y')
ylabel('z')
title('yz')
show()
