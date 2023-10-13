import numpy as np
import matplotlib.pyplot as plt
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


plot_density(0,0)
plot_density(1,0)
plot_density(0,0.5)
plot_density(1,0.5)