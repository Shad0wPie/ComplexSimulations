import lattice
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D  # NOQA


def visualize_simple_walk():
    l = lattice.Lattice(length=10,
                        num_total_individuals=1,
                        num_infected=1,
                        beta=0.5,
                        diffusion_rate=0.5,
                        gamma=0)
    l.run_simulation(plot_during=True, plot_final=False)


def detailed_plot_of_single_run():
    l = lattice.Lattice(length=100,
                        num_total_individuals=1000,
                        num_infected=10,
                        beta=0.5,
                        diffusion_rate=0.5,
                        gamma=0.005)
    l.run_simulation(plot_during=True, plot_final=True)


def _run_simulation(beta, r_0, diff_rate, iterations=5):
    gamma = np.float64(beta) / r_0
    print("running with beta: %s, r_0: %s, gamma: %s" % (beta, r_0, gamma))
    l = lattice.Lattice(length=100,
                        num_total_individuals=1000,
                        num_infected=10,
                        beta=beta,
                        diffusion_rate=diff_rate,
                        gamma=gamma)
    r_inf = 0
    for i in range(iterations):
        print("    iteration %s..." % i)
        l.reset()
        tmp_r_inf = l.run_simulation(plot_during=False, plot_final=False)
        print("        %s" % tmp_r_inf)
        r_inf += tmp_r_inf
    avg_r_inf = r_inf / iterations
    print("    average: %s" % avg_r_inf)
    return avg_r_inf


def plot_epidemic_thresholds(beta):
    diff_rate = 0.3
    r_0s = []
    r_infs = []
    for r_0 in np.linspace(0, 200, num=20):
        print("running with beta: %s and r_0: %s" % (beta, r_0))
        r_inf = _run_simulation(beta, r_0, diff_rate, iterations=5)
        print(r_inf)
        r_0s.append(r_0)
        r_infs.append(r_inf)
    print(r_0s, r_infs)
    plt.plot(r_0s, r_infs)
    plt.xlabel(r'$R_0$')
    plt.ylabel(r'$R_{\inf}$')
    plt.title(r'd=%s $\beta$=%s' % (diff_rate, beta))
    plt.show()


def plot_surface(beta_resolution, r0_resolution, iterations):
    beta_resolution *= 1j
    beta_resolution *= 1j
    BETA, R_0 = np.mgrid[0.01:1:beta_resolution, 0.01:250:r0_resolution]
    BETA = BETA.T
    R_0 = R_0.T

    print(R_0)
    print(BETA)

    vectorized_func = np.vectorize(_run_simulation, otypes=[np.float])

    R_INF = vectorized_func(BETA, R_0, 0.3, iterations=iterations)

    print(R_INF)

    # Plot the surface.
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(BETA, R_0, R_INF, cmap=cm.gray,
                           linewidth=0, antialiased=False)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$R_0$')
    ax.set_zlabel(r'$R_{\inf}$')
    ax.set_title(r'$\beta$-$R_0$-plane')

    plt.savefig('beta_r0_plane.pdf')
    plt.show()


if __name__ == '__main__':
    plot_epidemic_thresholds(beta=0.8)
    #    detailed_plot_of_single_run()
    # plot_surface(beta_resolution=10,r0_resolution=20, iterations=2)
