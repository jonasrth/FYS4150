import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

plt.rc("font", size=15)
plt.rc("xtick", labelsize=15)
plt.rc("ytick", labelsize=15)

def read_file(filename):
    """
    Returns columns in file as arrays
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

def read_cx_cube_file(filename):
    """
    Returns imaginary- and real part of data in text file as
    a 3D array
    """
    file = open("text_files/" + filename,"r")
    file.readline()
    s = file.readline().split()

    # checks if file is 2D array
    if len(s) == 2:
        M,N = int(s[0]),int(s[1])

        R = np.zeros((M,N))
        I = np.zeros_like(R)

        F = file.read().splitlines()

        for i in range(M):
            c = np.asarray([eval(s) for s in F[i].split()])
            R[i,:] = np.array(c[:,0])
            I[i,:] = np.array(c[:,1])
    else:
        M,N,S = int(s[0]),int(s[1]),int(s[2])

        R = np.zeros((S,M,N))
        I = np.zeros_like(R)

        F = file.read().splitlines()
        for n in range(S):
            for i in range(M):
                c = np.asarray([eval(s) for s in F[M*n+i].split()])
                R[n][i,:] = np.array(c[:,0])
                I[n][i,:] = np.array(c[:,1])

    file.close()

    return R,I

def plot_2DSE(data, time, filename, cbar_label=r"$p_{ij}$"):

    """
    Plots image of 2D Schrödinger equation results
    (wave function or probability distribution)
    -------
    data        - 2D array representing the data
    time        - time image was taken at
    filename    - name of file to save image to
    cbar_label  - label of colorbar
    """

    data = np.transpose(data)
    im = plt.imshow(data, extent=[0,1,0,1], origin="lower")
    cbar = plt.colorbar(im)
    cbar.set_label(cbar_label, fontsize=15)

    time_txt = plt.text(0.70, 0.92, "t = {:.3f}".format(time), color="white", fontsize=15)

    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("figures/"+filename)
    plt.close()



"""
Problem 7
"""

T = 0.008

# NO SLIT:

P = read_file("no_slit_total_prob.txt").flatten()
t = np.linspace(0,T,len(P))

plt.plot(t,abs(1-P))
plt.xlabel("t")
plt.ylabel(r"$|1 - p_{tot}(t)|$")
plt.savefig("figures/no_slit_total_prob.pdf")
plt.close()


# DOUBLE SLIT:

P = read_file("double_slit_total_prob.txt").flatten()
t = np.linspace(0,T,len(P))

plt.plot(t,abs(1-P))
plt.xlabel("t")
plt.ylabel(r"$|1 - p_{tot}(t)|$")
plt.savefig("figures/double_slit_total_prob.pdf")
plt.close()



"""
Problem 8
"""


R,I = read_cx_cube_file("double_slit_snapshots.txt")

P = R**2 + I**2

plot_2DSE(P[0],0, "p_double_slit_t0.pdf")
plot_2DSE(P[1],0.001, "p_double_slit_t1.pdf")
plot_2DSE(P[2],0.002, "p_double_slit_t2.pdf")

plot_2DSE(R[0],0, "Re_double_slit_t0.pdf", r"Re($u_{ij}$)")
plot_2DSE(R[1],0.001, "Re_double_slit_t1.pdf", r"Re($u_{ij}$)")
plot_2DSE(R[2],0.002, "Re_double_slit_t2.pdf", r"Re($u_{ij}$)")

plot_2DSE(I[0],0, "Im_double_slit_t0.pdf", r"Im($u_{ij}$)")
plot_2DSE(I[1],0.001, "Im_double_slit_t1.pdf", r"Im($u_{ij}$)")
plot_2DSE(I[2],0.002, "Im_double_slit_t2.pdf", r"Im($u_{ij}$)")



"""
Problem 9
"""


R,I = read_cx_cube_file("p_single.txt")

# Finding position of x = 0.8
i = int(np.round(0.8*R.shape[0]))
# Defining y-array
y = np.linspace(0,1,R.shape[0])


# SINGLE SLIT:

P = R**2 + I**2
plot_2DSE(P, 0.002, "p_single.pdf")

#plt.figure(figsize=(7,5))
plt.plot(y,P[i,:]/np.sum(P[i,:]))
plt.xlabel("y")
plt.ylabel(r"$p(x=0.8,y;t=0.002)$")
plt.savefig("figures/p_single_diff.pdf", bbox_inches="tight")
plt.close()


# DOUBLE SLIT:

R,I = read_cx_cube_file("p_double.txt")
P = R**2 + I**2
plot_2DSE(P, 0.002, "p_double.pdf")

#plt.figure(figsize=(8,5))
plt.plot(y,P[i,:]/np.sum(P[i,:]))
plt.xlabel("y")
plt.ylabel(r"$p(x=0.8,y;t=0.002)$")
plt.savefig("figures/p_double_diff.pdf", bbox_inches="tight")
plt.close()


# TRIPLE SLIT:

R,I = read_cx_cube_file("p_triple.txt")
P = R**2 + I**2
plot_2DSE(P, 0.002, "p_triple.pdf")

#plt.figure(figsize=(8,5))
plt.plot(y,P[i,:]/np.sum(P[i,:]))
plt.xlabel("y")
plt.ylabel(r"$p(x=0.8,y;t=0.002)$")
plt.savefig("figures/p_triple_diff.pdf", bbox_inches="tight")
plt.close()
