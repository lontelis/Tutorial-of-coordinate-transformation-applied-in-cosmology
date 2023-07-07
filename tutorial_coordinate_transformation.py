#############################################################################
#
# Tutorial of coordinate transformation applied in cosmology
# Author     : P.Ntelis 
# date       : 7 July 2023
# Comment 0:
# code found in repository: 
# https://github.com/lontelis/Tutorial-of-coordinate-transformation-applied-in-cosmology
# Comment 1
# dependences: cosmopit
# details found in: https://github.com/lontelis/cosmopit
# to get cosmopit just perform the following command in a terminal: 
# git clone https://github.com/lontelis/cosmopit.git
# Comment 2
# The tutorial is dedicated to a workshop during the conference of Cosmology from Home
# details found in: https://cosmologyfromhome.com/
#
#############################################################################

print('# import important modules')
import matplotlib.pyplot as plt
import numpy as np
import cosmopit
from cosmopit import galtools

print('# Fixing random state for reproducibility')
np.random.seed(19680801)


def randrange(n, vmin, vmax):
    """
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    """
    return (vmax - vmin)*np.random.rand(n) + vmin

n = 200
print('# producing random points in a 3D cartesian coordinate system')
xs = randrange(n, -1000, 1000)
ys = randrange(n, -1000, 1000)
zs = randrange(n, -1000, 1000)


print('# transform: x,y,z cartesian coordinates to radius, theta, phi spherical coordinates')
rs,ths,phs = galtools.xyz2rthph(xs,ys,zs)
print('# cosmology for homogeneous and isotropic universe')
print('# define cosmological parameters          ')
print(' Om  OL   w0  w1')
cosmopars=[0.3,0.7,-1.0,0.0]
print('# transform:  radius, theta, phi spherical coordinates to redshift (z), DEC, RA')
reds,decs,ras = np.zeros(len(xs)),np.zeros(len(xs)),np.zeros(len(xs))
for i in range(len(xs)):
	reds[i],decs[i],ras[i] = galtools.rthph2zdecra(rs[i],ths[i],phs[i],drho=0.0,dtheta=0.0,dphi=0.0,params=cosmopars)
print('# transform:  xyz cartesian coordinates to redshift (z), DEC, RA')
redsd,decsd,rasd = np.zeros(len(xs)),np.zeros(len(xs)),np.zeros(len(xs))
for i in range(len(xs)):
	redsd[i],decsd[i],rasd[i] = galtools.xyz2zdecra(xs[i],ys[i],zs[i],drho=0.0,dtheta=0.0,dphi=0.0,params=cosmopars)

print('# inverse transform: redshift (z), DEC, RA to radius, theta and phi spherical coordinates')
rsi,thsi,phsi = np.zeros(len(xs)),np.zeros(len(xs)),np.zeros(len(xs))
for i in range(len(xs)):
	rsi[i],thsi[i],phsi[i] = galtools.zdecra2rthph(reds[i],decs[i],ras[i],params=cosmopars,dist_type='proper')

print('# inverse transform:  radius, theta, phi spherical coordinates to cartesian coordinates x,y,z')
xsi, ysi, zsi = galtools.rthph2xyz(rs,ths,phs)

print('# inverse tansform: ')
xsid, ysid, zsid = np.zeros(len(xs)),np.zeros(len(xs)),np.zeros(len(xs))
for i in range(len(xs)):
	xsid[i], ysid[i], zsid[i] = galtools.zdecra2xyz(reds[i],decs[i],ras[i],params=cosmopars,dist_type='proper')

print('# tests of the routines of transformations:')
print("test 1) test_xyz2zdecra_zdecra2xyz")
galtools.test_xyz2zdecra_zdecra2xyz(x=xs, y=ys, z=zs, rtol=1e-4)
print("test 2) test_zdecra2rthph_rthph2zdecra")
galtools.test_zdecra2rthph_rthph2zdecra(z=reds,dec=decs,ra=ras,rtol=1e-3)
print("test 3) test_xyz2rthph_rthph2xyz")
galtools.test_xyz2rthph_rthph2xyz(x=xs, y=ys, z=zs,rtol=1e-3)
print("test 4) test_rthph2xyz_xyz2rthph")
for i in range(len(xs)):
	galtools.test_rthph2xyz_xyz2rthph(r=rs[i],th=ths[i],ph=phs[i])

print('# Plotting functions')


plt.ion()
fig = plt.figure(1,figsize=(10,6))
plt.suptitle('Cartesian coordinates 3D: $(x,y,z)$')
ax = fig.add_subplot(projection='3d')
ax.scatter(xs, ys, zs, marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show(),plt.ion()


fig = plt.figure(2,figsize=(10,6))
plt.suptitle('Cartesian coordinates 3D: $(x,y,z)$')
ax = fig.add_subplot(projection='3d')
ax.scatter(rs, ths, phs, marker='o')
ax.set_xlabel('$r$ Label')
ax.set_ylabel('$\\theta$ Label')
ax.set_zlabel('$\\phi$ Label')
plt.show(),plt.ion()


fig = plt.figure(3,figsize=(10,6))
plt.suptitle('Astronomical coordinates 3D: (redshift, Dec, R.A.) ')
ax = fig.add_subplot(projection='3d')
ax.scatter(reds, decs, ras, marker='o')
ax.set_xlabel('$z$ Label')
ax.set_ylabel('$DEC$ Label')
ax.set_zlabel('$RA$ Label')
plt.show(),plt.ion()

plt.figure(4,figsize=(10,6))
plt.clf()
plt.suptitle('Cartesian coordinates 3D, 2D profiles: $(x,y,z)$')
plt.subplot(121)
plt.plot(xs,ys,'o')
plt.grid()
plt.xlabel('X label')
plt.ylabel('Y label')
plt.subplot(122)
plt.plot(xs,zs,'o')
plt.grid()
plt.xlabel('X label')
plt.ylabel('Z label')
plt.show()


plt.figure(5,figsize=(10,6))
plt.clf()
plt.suptitle('Spherical coordinates 3D, 2D profiles: $(r,\\theta,\\phi)$')
plt.subplot(121)
plt.plot(rs,ths,'o')
plt.grid()
plt.xlabel('r label')
plt.ylabel('$\\theta$ label')
plt.subplot(122)
plt.plot(rs,phs,'o')
plt.grid()
plt.xlabel('r label')
plt.ylabel('$\\phi$ label')
plt.show()

plt.figure(6,figsize=(10,6))
plt.clf()
plt.suptitle('Astronomical coordinates 3D, 2D profiles: (redshift, Dec, R.A.) ')
plt.subplot(221)
plt.plot(reds,ras,'o')
plt.grid()
plt.xlabel('$z$ label')
plt.ylabel('$R.A.$ label')
plt.subplot(222)
plt.plot(ras,decs,'o')
plt.grid()
plt.xlabel('$R.A.$ label')
plt.ylabel('$Dec$ label')
plt.subplot(223)
plt.plot(reds,decs,'o')
plt.grid()
plt.xlabel('$z$ label')
plt.ylabel('$Dec$ label')
plt.show()


plt.figure(7,figsize=(10,6))
plt.clf()
plt.title('radial comoving distance vs redshift')
plt.plot(reds,rs,'.')
plt.grid()
plt.xlabel('$z$, redshift')
plt.ylabel('$r_C(z)$, radial comoving distance [Mpc/h]')
plt.show()


plt.figure(8,figsize=(10,6))
plt.clf()
plt.plot(reds,reds,'k--')
plt.plot(reds,redsd,'.')
plt.grid()
plt.xlabel('$z$, redshift from two functions writen separatly in this code')
plt.ylabel('$z$, redshift from two function built in the same module')
plt.show()

