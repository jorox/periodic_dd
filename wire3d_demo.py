'''
=================
3D wireframe plot
=================

A very basic demonstration of a wireframe plot.
'''

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

def u_dip(X,Y,b,R):
    """
    return the displacement field of a dislocation in the range x,y
    The dipole dislocations are at y=0 x = R[0] with Burgers b
    and x=R[1] with burgers -b
    """
    if X.shape != Y.shape: 
        print("ERROR")
        return 0
    Z = np.zeros(X.shape)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Z[i,j] = b/2.0/3.1415*(np.arctan2(Y[i,j]-R[0][1],X[i,j]-R[0][0])
                    -np.arctan2(Y[i,j]-R[1][1],X[i,j]-R[1][0]))
    return Z


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 1. Define the full region (zero-th image + periodic)
x_range = np.arange(-10,10.1,0.3); N = len(x_range)
y_range = np.arange(-5,5.1,0.03); M = len(y_range)
# 2. Define the region of interest i.e. the zero-th image
roi = [[-1,1],[-0.5,0.5]]

# 3. Define the 2D matrices for calculation of u_z
X = []
Y = []
for i in range(M): X.append(x_range) # X = MxN
for i in range(N): Y.append(y_range) # Y = NxM
X = np.array(X)
Y = np.transpose(np.array(Y)) # Y = MxN

# u_dipole
Z = u_dip(X,Y,1,[[0.5,0.0],[-0.5,0.0]])

# Find centers of dipoles in periodic images
R_0 = np.arange(x_range[0],x_range[-1],roi[0][1]-roi[0][0])
R_1 = np.arange(y_range[0],y_range[-1],roi[1][1]-roi[1][0])
print(R_0)
print(R_1)
R = []
for i in range(len(R_0)):
    for j in range(len(R_1)):
        if R_0[i]==0 and R_1[j]==0: continue
        R.append([R_0[i],R_1[j]])

# Calculate u_sum
for i in range(len(R)):
    Z += u_dip(X,Y,1,[[R[i][0]+0.5,R[i][1]], [R[i][0]-0.5,R[i][1]]])

# Calculate u_err and u_pbc
d = np.zeros((3,3))


for i in range(M):
    for j in range(N):
        if (X[i,j]>roi[0][1] or X[i,j]<roi[0][0] or Y[i,j]>roi[1][1] or Y[i,j]<roi[1][0]):
            Z[i,j] = np.NaN
        else:
            pass


# Plot a basic wireframe.
ax.set_xlim3d([-1 , 1])
ax.set_ylim3d([-0.5,0.5])
ax.set_zlim3d([-1, 1])
ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

plt.show()
