import numpy as np
import matplotlib.pyplot as plt

# Makes figures smaller, but makes text larger in comparison (easier to read):
plt.rc('figure', figsize=(5.5, 4))


def read_file(filename):
    """
    Returns arrays of columns in file
    """
    file = open("text_files/" + filename,"r")
    file.readline()
    file.readline()

    F = file.read().splitlines()
    for i in range(len(F)):
        F[i] = F[i].split()
    F = np.array(F)
    F = F.astype(float)

    file.close()

    return F


F = read_file("2p_positions.txt")

t  = F[:,0]
R1 = F[:,1:4]
V1 = F[:,4:7]
R2 = F[:,7:10]
V2 = F[:,10:13]


# Plotting positions in xy plane for two particles without interactions:

plt.plot(R1[:,0],R1[:,1],label="particle 1")
plt.plot(R2[:,0],R2[:,1],label="particle 2")
plt.legend()
plt.axis("equal")
plt.xlabel(r"x [$\mu$m]")
plt.ylabel(r"y [$\mu$m]")
plt.savefig("figures/2p_xy.pdf")
plt.show()


# Plotting z position over time

plt.plot(t,R1[:,2])
plt.axis("equal")
plt.xlabel(r"t [$\mu$s]")
plt.ylabel(r"z [$\mu$m]")
plt.savefig("figures/1p_tz.pdf")
plt.show()


# Phase plots without interactions

fig,ax = plt.subplots(1,3)
plots = [("x","r"),("y","b"),("z","g")]
for i in range(3):
    ax[i].plot(R1[:,i],V1[:,i],color=plots[i][1])
    ax[i].set_xlabel(rf"${plots[i][0]}$")
    ax[i].set_ylabel(rf"$v_{plots[i][0]}$")
plt.tight_layout()
plt.savefig("figures/phase_plot.pdf")
plt.show()


# Relative error plots for Forward Euler

for i in range(5):
    t_rel = read_file("FE_rel_error_h"+str(i)+".txt")
    plt.plot(t_rel[1:,0],t_rel[1:,1],label=r"$h = 10^{-%s}\mu$s" % i)
plt.legend()
plt.yscale("log")
plt.xlabel("t [$\mu$s]")
plt.ylabel("relative error")
plt.savefig("figures/FE_rel_error.pdf")
plt.show()

# Relative error plots for Runge-Kutta 4

for i in range(5):
    t_rel = read_file("RK4_rel_error_h"+str(i)+".txt")
    plt.plot(t_rel[1:,0],t_rel[1:,1],label=r"$h = 10^{-%s}\mu$s" % i)
plt.legend()
plt.yscale("log")
plt.xlabel("t [$\mu$s]")
plt.ylabel("relative error")
plt.savefig("figures/RK4_rel_error.pdf")
plt.show()


F = read_file("2p_interaction_positions.txt")

t  = F[:,0]
R1 = F[:,1:4]
V1 = F[:,4:7]
R2 = F[:,7:10]
V2 = F[:,10:13]

# xy plane plot for two particles, with interactions

plt.plot(R1[:,0],R1[:,1],label="particle 1")
plt.plot(R2[:,0],R2[:,1],label="particle 2")
plt.legend()
plt.axis("equal")
plt.xlabel(r"x [$\mu$m]")
plt.ylabel(r"y [$\mu$m]")
plt.savefig("figures/2p_interaction_xy.pdf")
plt.show()

# Phase plots with interactions

fig,ax = plt.subplots(1,3)
plots = [("x","r"),("y","b"),("z","g")]
for i in range(3):
    ax[i].plot(R1[:,i],V1[:,i],color=plots[i][1])
    ax[i].set_xlabel(rf"${plots[i][0]}$")
    ax[i].set_ylabel(rf"$v_{plots[i][0]}$")
plt.tight_layout()
plt.savefig("figures/interaction_phase_plot.pdf")
plt.show()



"""
Particles lost as function of induced frequency
"""


NW = read_file("NW_f0.1.txt")
omegaV = NW[:,0]
n = NW[:,1]
plt.plot(omegaV,n,label="$f = 0.1$")

NW = read_file("NW_f0.4.txt")
omegaV = NW[:,0]
n = NW[:,1]
plt.plot(omegaV,n,label="$f = 0.4$")


NW = read_file("NW_f0.7.txt")
omegaV = NW[:,0]
n = NW[:,1]
plt.plot(omegaV,n,label="$f = 0.7$")

plt.legend()
plt.xlabel(r"$\omega_V [\mu$s$^{-1}]$")
plt.ylabel("escaped particles / total particles")
plt.savefig("figures/NW.pdf")
plt.show()

"""
zooming in on peak:
"""

NW = read_file("NW_f0.1_zoom.txt")

omegaV = NW[:,0]
n = NW[:,1]
plt.plot(omegaV,n,label="interactions off")

NW = read_file("NW_f0.1_zoom_interaction.txt")

omegaV = NW[:,0]
n = NW[:,1]
plt.plot(omegaV,n,label="interactions on")

plt.legend()
plt.xlabel(r"$\omega_V [\mu$s$^{-1}]$")
plt.ylabel("escaped particles / total particles")
plt.savefig("figures/NW_zoom.pdf")
plt.show()
