import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

plt.figure(figsize=(12,4))
ap = ArgumentParser()
ap.add_argument("--filename","-f",type=str,help="basis file name")
ap.add_argument("--epsilon", "-e", type=float,default=0.5, help="fractional flux value")

aa = ap.parse_args()

words = {1: "one_pin", 2: "two_pins"}

with open(f"berry_phase/{aa.filename}_two_pins_ϵ_{aa.epsilon}.txt") as f:
	data = [x.split() for x in f.readlines()]

R = np.array([float(x[0]) for x in data])
D = [float(x[1]) for x in data]
P = np.array([float(x[2])/(2*3.14159265) for x in data])

with open(f"berry_phase/{aa.filename}_two_pins_Lz_ϵ_{aa.epsilon}.txt") as f:
	data = [x.split() for x in f.readlines()]

L = np.array([-np.mod(float(x[1]),1) for x in data])


theory = np.mod(-(2*R)**2,1)-1 # -πR²/2π × 2
print("Difference between phase and area:")
br = P-theory
print(br)

plt.plot(R,br, "o",markerfacecolor="none", label="double exchange phase")

plt.xlabel(r"$R$")
plt.ylabel(r"Berry phase / 2$\pi$")
plt.legend()
plt.savefig(f"{aa.filename}_two_pins_ϵ_{aa.epsilon}.svg")