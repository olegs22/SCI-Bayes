import numpy as np
import matplotlib.pyplot as plt
import datetime
from sidereal import *
from times import *
import seaborn as sb
import healpy as hp
sb.set_style('white')
rad = np.pi/180.

direct = '/Users/Oleg/Documents/SCI-Bayes/avances/'
data = np.loadtxt(direct+'70MHz.txt')

glon = -118.3011 * rad
glat = float("{0:.4f}".format(28.9733 * rad))
theta = data[:,0]*rad
phi = data[:,1]*rad
temp = data[:,2]

phi_range = len(np.arange(0,91,1))
#nblocks = len(phi) / phi_range #this is the number of blocks i'll make

n_theta=[]
n_phi=[]
n_data=[]
for i in range(len(theta)):
    if theta[i] > 0.0:
        n_theta.append(theta[i])
        n_phi.append(phi[i])
        n_data.append(temp[i])

nblocks = len(n_phi) / phi_range #this is the number of blocks i'll make
print nblocks

def rotations(theta,phi,psi):
    rot_z = hp.rotator.Rotator(rot=(phi, 0, 0), deg=True, eulertype='y')
    rot_y = hp.rotator.Rotator(rot=(0, -theta, 0), deg=True, eulertype='y')
    rot_zp = hp.rotator.Rotator(rot=(psi,0,0),deg=True,eulertype='y')
    return rot_z * rot_y * rot_zp


rotator = rotations(glat/rad,glon/rad,0.0)
rotator = rotator.get_inverse()

#nside = 512/4
#bmap = np.zeros(hp.nside2npix(nside), dtype=np.double)
#grid = np.ndarray((trans[0,:],trans[1,:]))
#len(trans[0,:])
#len(trans[1,:])
rot_theta = []
rot_phi = []
for i in range(len(n_theta)):
    for k in range(len(n_phi)):
        (rtheta,rphi) = rotator((n_theta[i]/rad,n_phi[k]/rad))
    rot_theta.append(rtheta*rad)
    rot_phi.append(rphi*rad)
"""
dec = rot_theta

ra = []
for i in range(len(rot_phi)):
    val = datetime.datetime(2013,6,14,13,59,56)
    valt = SiderealTime.fromDatetime(val)
    LSTT = SiderealTime.lst(valt,rot_phi[i])*15.*rad
    ra.append(LSTT)


alpha_g=192.85*rad
delta_g=27.128333*rad
alpha_c=266.4*rad
delta_c=-28.929656*rad

trans = []
for i in range(len(ra)):
    B=np.arcsin(np.sin(dec[i])*np.sin(delta_g)+np.cos(dec[i])*np.cos(delta_g)*np.cos(ra[i]-alpha_g))
    J=(np.sin(dec[i])*np.cos(delta_g)-np.cos(dec[i])*np.sin(delta_g)*np.cos(ra[i]-alpha_g))/np.cos(B)
    K=np.arcsin(np.cos(dec[i])*np.sin(ra[i]-alpha_g)/np.cos(B))
    Q=np.arccos(np.sin(delta_c)/np.cos(delta_g))

    if J<0.:
        L=Q/rad+K/rad-180.
    else:
        L=(Q/rad-K/rad)

    if L<0.:
        L=L+360.

    trans.append(np.array([L*rad, B]))
trans = np.transpose(trans)


for i in range(len(trans[0,:])):
    if trans[0,i] > np.pi:
        trans[0,i] -= np.pi

nside = 512/4
pix = hp.ang2pix(nside, trans[0,:], trans[1,:])
#pix2 = hp.ang2pix(nside, ra, dec)
#pix = hp.ang2pix(nside, rotated_grid[0], rotated_grid[1])
bmap = np.zeros(hp.nside2npix(nside), dtype=np.double)
#bmap2 = np.zeros(hp.nside2npix(nside), dtype=np.double)
bmap[pix] = n_data
#bmap2[pix2] = data[:,2]
hp.mollview(np.log10(bmap))
plt.show()
#gx_data = np.loadtxt('70MHz.dat')

#covolution = np.sum(bmap * gx_data) / np.sum(bmap)
#print covolution
"""
