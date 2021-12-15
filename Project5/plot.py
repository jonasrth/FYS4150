import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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

def read_cube_file(filename):
    """
    Returns columns in file as arrays
    """
    file = open("text_files/" + filename,"r")
    file.readline()
    s = file.readline().split()
    M,N,S = int(s[0]),int(s[1]),int(s[2])

    C = np.zeros((S,M,N))


    F = file.read().splitlines()
    for n in range(S):
        for i in range(M):
            C[n][i,:] = np.asarray(F[M*n+i].split(),dtype="float")

    file.close()

    return C

def read_cx_cube_file(filename):
    """
    Returns columns in file as arrays
    """
    file = open("text_files/" + filename,"r")
    file.readline()
    s = file.readline().split()
    M,N,S = int(s[0]),int(s[1]),int(s[2])

    R = np.zeros((S,M,N))
    I = np.zeros_like(R)

    R_ = np.zeros((M,N))
    I_ = np.zeros_like(R_)


    F = file.read().splitlines()
    for n in range(S):
        for i in range(M):
            c = np.asarray([eval(s) for s in F[M*n+i].split()])
            R[n][i,:] = np.array(c[:,0])
            I[n][i,:] = np.array(c[:,1])

    file.close()

    return R,I

"""
F = np.transpose(read_file("test1.txt"))
#F = read_file("test.txt")

im = plt.imshow(F,origin="lower")
plt.colorbar(im)

plt.show()

F = np.transpose(read_file("test2.txt"))
#F = read_file("test.txt")

im = plt.imshow(F,origin="lower")
plt.colorbar(im)

plt.show()

F = np.transpose(read_file("test3.txt"))
#F = read_file("test.txt")

im = plt.imshow(F,origin="lower")
plt.colorbar(im)

plt.show()
"""


def Animate2DSE(z_data_list,dt,T,file=None):

    t_points = np.linspace(0,T,z_data_list.shape[0])

    #
    # Now the list z_data_list contains a series of "frames" of z(x,y,t),
    # where each frame can be plotted as a 2D image using imshow. Let's
    # animate it!
    #

    # Some settings
    fontsize = 12
    t_min = t_points[0]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

    # Plot the first frame
    img = ax.imshow(np.transpose(z_data_list[0]), extent=[0,1,0,1], cmap=plt.get_cmap("viridis"), norm=norm, origin="lower")

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(r"$p(x,y;t)$", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.65, 0.95, "t = {:.3e}".format(t_min), color="white", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
        img.set_norm(norm)

        # Update z data
        img.set_data(np.transpose(z_data_list[i]))

        # Update the time label
        current_time = t_min + i*dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

    # Run the animation!
    plt.show()

    ## Save the animation
    if file != None:
        anim.save("./"+file+".mp4", writer="ffmpeg", bitrate=-1, fps=30)

#R,I = read_cx_cube_file("cube.txt")
#z_data_list = R**2 + I**2

z_data_list = read_cube_file("cube.txt")

Animate2DSE(z_data_list,2.5e-5,0.002,file="test2")


"""
fig,ax = plt.subplots()

ax.plot(abs(1-np.sum(z_data_list,axis=(1,2))))
ax.set_yscale("log")
plt.show()

fig,ax = plt.subplots()

M = len(z_data_list[0,0,:])
m = int(round(0.8*M))

p = z_data_list[-1,m,:]

ax.plot(p/np.sum(p))
plt.show()
"""
