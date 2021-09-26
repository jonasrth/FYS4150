import numpy as np
import matplotlib.pyplot as plt

"""
N = np.arange(2,21)
iterations = np.asarray([1,10,6,37,46,80,108,140,179,228,275,315,348,424,470,550,619,703,767])

g = np.gradient(np.log(iterations),np.log(N))
#g = np.gradient(np.log(N),np.log((767/400)*N**2))

#plt.plot(N,iterations)
#plt.plot(N,(767/400)*N**2)
#plt.plot(np.log(N),np.log(iterations))
plt.plot(N,g)
plt.show()
"""

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

legend_elements_v = [r"$v_{%s}(\hat x)$" % i for i in range(3)]
legend_elements_u = [r"$u_{%s}(\hat x)$" % i for i in range(3)]
legend_elements = legend_elements_v+legend_elements_u

# Problem 7a:

x,v = read_file("solution_N10.txt")
x,u = read_file("solution_N10_analytical.txt")

fig,ax = plt.subplots()
ax.plot(x,v)
ax.plot(x,u)
ax.set_xlabel(r"$\hat x$")
ax.legend(legend_elements)
ax.set_title("Comparing analytical- and numerical solutions for $N = 10$")
plt.savefig("figures/comparison_N10.pdf")
plt.show()


# Problem 7b:

x,v = read_file("solution_N100.txt")
x,u = read_file("solution_N100_analytical.txt")

fig,ax = plt.subplots()
ax.plot(x,v)
ax.plot(x,u)
ax.set_xlabel(r"$\hat x$")
ax.legend(legend_elements)
ax.set_title("Comparing analytical- and numerical solutions for $N = 100$")
plt.savefig("figures/comparison_N100.pdf")
plt.show()
