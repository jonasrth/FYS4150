import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.stats import linregress

plt.rcParams['font.size'] = '14'

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

# Analytical partition function, expectation values, heat capacity
# and magnetic susceptibility for LxL lattice with length L = 2 as
# a function of temperature T.
Z_L2 = lambda T: 2*np.exp(8/T) + 12 + 2*np.exp(-8/T)
eps_L2 = lambda T: 4/Z_L2(T)*(np.exp(-8/T) - np.exp(8/T))
eps2_L2 = lambda T: 8/Z_L2(T)*(np.exp(-8/T) + np.exp(8/T))
m_L2 = lambda T: 2/Z_L2(T)*(np.exp(8/T) + 2)
m2_L2 = lambda T: 2/Z_L2(T)*(np.exp(8/T) + 1)
CV_L2 = lambda T: 4/T**2*(eps2_L2(T) - eps_L2(T)**2)
chi_L2 = lambda T: 4/T*(m2_L2(T) - m_L2(T)**2)

print("<eps>   [J]   : " + str(eps_L2(1.0)))
print("<|m|>   [1]   : " + str(m_L2(1.0)))
print("<eps^2> [J^2] : " + str(eps2_L2(1.0)))
print("<m^2>   [1]   : " + str(m2_L2(1.0)))
print("C_V     [k_B] : " + str(CV_L2(1.0)))
print("chi     [1/J] : " + str(chi_L2(1.0)))

fig1,ax1 = plt.subplots()
fig2,ax2 = plt.subplots()

n = np.arange(len(read_file("EM_T1_ordered_L20.txt")[:,0]))+1

EM = np.cumsum(read_file("EM_T1_unordered_L20.txt"),axis=0)/n[:,np.newaxis]
ax1.plot(n,EM[:,0])
ax2.plot(n,EM[:,1])

EM = np.cumsum(read_file("EM_T1_ordered_L20.txt"),axis=0)/n[:,np.newaxis]
ax1.plot(n,EM[:,0])
ax2.plot(n,EM[:,1])

EM = np.cumsum(read_file("EM_T2_unordered_L20.txt"),axis=0)/n[:,np.newaxis]
ax1.plot(n,EM[:,0])
ax2.plot(n,EM[:,1])

EM = np.cumsum(read_file("EM_T2_ordered_L20.txt"),axis=0)/n[:,np.newaxis]
ax1.plot(n,EM[:,0])
ax2.plot(n,EM[:,1])

ax1.legend([r"$T = 1.0$ $J/k_B$, unordered",r"$T = 1.0$ $J/k_B$, ordered",
            r"$T = 2.4$ $J/k_B$, unordered",r"$T = 2.4$ $J/k_B$, ordered"])
ax2.legend([r"$T = 1.0$ $J/k_B$, unordered",r"$T = 1.0$ $J/k_B$, ordered",
            r"$T = 2.4$ $J/k_B$, unordered",r"$T = 2.4$ $J/k_B$, ordered"])
ax1.set_xlabel("MC cycles")
ax2.set_xlabel("MC cycles")
ax1.set_ylabel(r"$<\epsilon>$")
ax2.set_ylabel(r"$<|m|>$")
fig1.savefig("figures/eps_L20.pdf")
fig2.savefig("figures/m_L20.pdf")
plt.show()


# burn-in time:
n_b = 2000
step = 4/20/20

E = read_file("EM_T1_L20.txt")[n_b:,0]

bins = int((np.max(E)-np.min(E))/step)
hist, bin_edges = np.histogram(E,bins = bins)
weights = np.ones_like(E)/len(E)

fig,ax = plt.subplots()
ax.hist(E,bins=bins,weights=weights,align="left",rwidth=0.85)
ax.set_xlabel(r"$\epsilon$ [J]")
ax.set_ylabel(r"$p(\mathbf{s}$; $T = 1.0$ $J/k_B)$")
fig.savefig("figures/EM_T1_L20.pdf")
plt.show()



E = read_file("EM_T2_L20.txt")[n_b:,0]

bins = int((np.max(E)-np.min(E))/step)
hist, bin_edges = np.histogram(E,bins = bins)

weights = np.ones_like(E)/len(E)
fig,ax = plt.subplots(figsize = (8,5))
ax.hist(E,bins=bins,weights=weights,align="left",rwidth=0.85)
ax.set_xlabel(r"$\epsilon$ [J]")
ax.set_ylabel(r"$p(\mathbf{s}$; $T = 2.4$ $J/k_B)$")
fig.savefig("figures/EM_T2_L20.pdf")
plt.show()


L_values = [40,60,80,100]

dict = {}


for L in L_values:
    dict[L] = read_file("init_temps_L"+str(L)+".txt")

for L in L_values:
    T = dict[L][:,0]
    eps = dict[L][:,1]
    plt.plot(T,eps,label="L = %d" % L)
plt.legend()
plt.xlabel(r"$T$ [J/k$_B$]")
plt.ylabel(r"$<\!\varepsilon\!>$ [J]")
plt.savefig("figures/eps_T.pdf")
plt.show()

for L in L_values:
    T = dict[L][:,0]
    m = dict[L][:,2]
    plt.plot(T,m,label="L = %d" % L)
plt.legend()
plt.xlabel(r"$T$ [J/k$_B$]")
plt.ylabel(r"$<\!|m|\!>$ [1]")
plt.savefig("figures/m_T.pdf")
plt.show()

for L in L_values:
    T = dict[L][:,0]
    C_V = dict[L][:,3]
    plt.plot(T,C_V,label="L = %d" % L)
plt.legend()
plt.xlabel(r"$T$ [J/k$_B$]")
plt.ylabel(r"$C_V$ [k$_B$]")
plt.savefig("figures/CV_T.pdf")
plt.show()

for L in L_values:
    T = dict[L][:,0]
    chi = dict[L][:,4]
    plt.plot(T,chi,label="L = %d" % L)
plt.legend()
plt.xlabel(r"$T$ [J/k$_B$]")
plt.ylabel(r"$\chi$ [1]")
plt.savefig("figures/chi_T.pdf")
plt.show()


for L in L_values:
    dict[L] = read_file("zoom_peak_temps_2_L"+str(L)+".txt")

T_c = np.zeros(len(L_values))
i = 0
for L in L_values:
    T = dict[L][:,0]
    C_V = dict[L][:,3]
    C_V_ = savgol_filter(C_V,41,3)
    T_c[i] = T[np.argmax(C_V_)]
    plt.plot(T,C_V,label="L = %d" % L)
    plt.plot(T,C_V_,"C"+str(i)+"--")
    i += 1
plt.legend()
plt.xlabel(r"$T$ [J/k$_B$]")
plt.ylabel(r"$C_V$ [k$_B$]")
plt.savefig("figures/CV_T2.pdf")
plt.show()

fig,ax = plt.subplots(figsize=(8,5))
ax.plot(L_values,T_c)
ax.set_xlabel("L")
ax.set_ylabel(r"$T_c$ [$J/k_B$]")
fig.savefig("figures/T_c.pdf")
plt.show()

L = np.asarray(L_values)
fig,ax = plt.subplots()
ax.plot(L,T_c*L,label=r"$T_c(L=\infty)L + a$")
ax.set_xlabel("L")
ax.set_ylabel(r"$T_cL$ [$J/k_B$]")
ax.legend()
fig.savefig("figures/T_cL.pdf")
plt.show()

res = linregress(L,T_c*L)
print("T_C(L = infinity) = " + str(res.slope)+ " J/k_B")
