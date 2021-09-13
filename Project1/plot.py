import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    """
    Returns arrays of columns in file
    """
    file = open("text_files/"+filename,"r")
    file.readline()

    x = np.array([])
    y = np.array([])

    for line in file:
        w = line.split()

        x = np.append(x,float(w[0]))
        y = np.append(y,float(w[1]))

    return x,y

# Problem 2:

x,u = read_file("u_exact.txt")
plt.plot(x,u)
plt.title(r"Plot of exact solution to $u(x)$")
plt.xlabel(r"$x$")
plt.ylabel(r"$u(x)$")
plt.savefig("figures/u_exact.pdf")
plt.show()

# Problem 7:

plt.plot(x,u,label="u(x)")
for n in [10,100,1000]:
    x,v = read_file("v_n"+str(n)+".txt")
    plt.plot(x,v,label="n = "+str(n))
x,u = read_file("u_exact.txt")
plt.legend()
plt.title(r"Plot of approximated $u(x)$ for different $n$'s next to exact solution")
plt.xlabel(r"$x$")
plt.ylabel(r"$v(x)$")
plt.savefig("figures/u_approximate.pdf")
plt.show()

# Problem 8a:

for n in [10,100,1000,10000,100000]:
    x,abs_error = read_file("abs_error_n"+str(n)+".txt")
    plt.plot(x,np.log10(abs_error),label="n = "+str(n))
plt.legend()
plt.title(r"Logarithm of the absolute error as a function of $x$, for different $n$'s")
plt.xlabel(r"$x$")
plt.ylabel(r"log$_{10}(\Delta_i)$")
plt.savefig("figures/abs_error.pdf")
plt.show()

# Problem 8b:

for n in [10,100,1000,10000,100000]:
    x,rel_error = read_file("rel_error_n"+str(n)+".txt")
    plt.plot(x,np.log10(rel_error),label="n = "+str(n))
plt.legend()
plt.title(r"Logarithm of the relative error as a function of $x$, for different $n$'s")
plt.xlabel(r"$x$")
plt.ylabel(r"log$_{10}(\epsilon_i)$")
plt.savefig("figures/rel_error.pdf")
plt.show()

# Problem 8c:

x,y = read_file("max_error.txt")
plt.plot(np.log10(x),np.log10(y))
plt.title(r"Maximum relative error as a function of $n$")
plt.xlabel(r"log$_{10}(n)$")
plt.ylabel(r"log$_{10}[$max$(\epsilon_i)]$")
plt.savefig("figures/max_error.pdf")
plt.show()

# Problem
