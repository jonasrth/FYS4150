import numpy as np
import matplotlib.pyplot as plt

set_matplotlib_formats('svg')
plt.rc('legend', frameon=False)
plt.rc('figure', figsize=(7, 7 / 1.75)) # Larger figure sizes
plt.rc('font', size=12)

def read_file(filename):
    """
    Returns arrays of columns in file
    """
    file = open("text_files/"+filename,"r")
    file.readline()

    F = file.read().splitlines()
    for i in range(len(F)):
        F[i] = F[i].split()
    F = np.array(F)
    F = F.astype(float)

    file.close()

    x = F[:,0]
    y = F[:,1:]

    return x,y

# Problem 6a:

N,iters = read_file("transformations_tridiag.txt")
fig,ax = plt.subplots()
ax.plot(N,iters)
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel("Transformations")
ax.set_title(r"Transformations for $N\times N$ tridiagonal matrix")
plt.savefig("figures/transformations_tridiag.pdf")
plt.show()

fig,ax = plt.subplots()
ax.plot(N,np.gradient(np.log(iters[:,0]),np.log(N)))
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\alpha$")
ax.set_title(r"Gradient of the logarithm of the number of transformations, over $N$")
plt.savefig("figures/gradient_tridiag.pdf")
plt.show()

# Problem 6b:

fig,ax = plt.subplots()
ax.plot(N,iters,label="tridiagonal")
N,iters = read_file("transformations_dense.txt")
ax.plot(N,iters,label="dense")
ax.legend()
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel("Transformations")
ax.set_title(r"Transformations for $N\times N$ dense matrix")
plt.savefig("figures/transformations_dense.pdf")
plt.show()


# Problem 7:

legend_elements_v = [r"$v_{%s}(\hat x)$" % i for i in range(3)]
legend_elements_u = [r"$u_{%s}(\hat x)$" % i for i in range(3)]
legend_elements = legend_elements_v+legend_elements_u

# Problem 7a:

x,v = read_file("solution_N10.txt")
x,u = read_file("solution_N10_analytical.txt")

fig,ax = plt.subplots()
ax.plot(x,v)
ax.plot(x,u)
ax.legend(legend_elements)
ax.grid()
ax.set_xlabel(r"$\hat x$")
ax.set_title("Comparing analytical- and numerical solutions for $N = 10$")
plt.savefig("figures/comparison_N10.pdf")
plt.show()


# Problem 7b:

x,v = read_file("solution_N100.txt")
x,u = read_file("solution_N100_analytical.txt")

fig,ax = plt.subplots()
ax.plot(x,v)
ax.plot(x,u)
ax.legend(legend_elements)
ax.grid()
ax.set_xlabel(r"$\hat x$")
ax.set_title("Comparing analytical- and numerical solutions for $N = 100$")
plt.savefig("figures/comparison_N100.pdf")
plt.show()
