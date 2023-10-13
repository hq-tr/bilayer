import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from time import time

def single_particle_function(m,s):
	return lambda z: z**(m+s) * np.exp(-abs(z)**2)

def density_distribution(positions, m, s):
	wf = single_particle_function(m,s)
	return abs(wf(positions))**2

def plot_density(m,s):
	st = time()
	x = np.linspace(-5,5,100)
	y = np.linspace(-5,5,100)
	X, Y = np.meshgrid(x,y)
	U = X + 1j * Y 
	D = density_distribution(U,m,s)
	plt.imshow(D)
	plt.savefig(f"density_m_{m}_shift_{s}.svg")
	plt.title(f"m = {m}, shift = {s}")
	print(f"{time()-st} seconds.")
	return

sqfactorial = lambda x: np.prod(np.sqrt(np.arange(1,x+1)))

def translation_coefficient(a,m):
	return np.prod(a/np.arange(1,m+1)) # The factor 1/10 is to avoid big number in the factorial
	
# Params:
No = 16
xp  = 3
yp  = 0
a  = (xp + 1j*yp)

s = 0.5 # Shift: either 0 or 0.5

# Calculation tools
coefs = np.array([translation_coefficient(a, m) for m in range(No+1)])
coefs /= np.sqrt(np.abs(np.dot(np.conj(a),a)))
wfs   = [single_particle_function(m, s) for m in range(No+1)]

# Prepare plot surface
Rmax = np.sqrt(2*No)
N = 100 # Number of points along each axis
x = np.linspace(-Rmax,Rmax,N)
y = np.linspace(-Rmax,Rmax,N)
X, Y = np.meshgrid(x,y)
U = (X + 1j * Y)/2
D = np.zeros((N,N))

st = time() # Time the for loop
for i in range(No):
	D += abs(coefs[i])**2 * abs(wfs[i](U))**2 # diagonal terms
	for j in range(i):
		term = np.conj(wfs[i](U)) * wfs[j](U) * np.conj(coefs[i]) * coefs[j]
		D += 2*np.real(term)

print(f"{time()-st} seconds.")
plt.pcolormesh(X,Y,D,edgecolor="none")
plt.axis("equal")
plt.plot([3],[0],"ro")
plt.title(f"x = {xp:.2f}, y = {yp:.2f}, shift = {s}")
plt.savefig(f"density_x_{xp:.2f}_y_{yp:.2f}_shift_{s}.svg")

