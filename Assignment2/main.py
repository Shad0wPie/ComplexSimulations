import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.measurements import label
from scipy.optimize import curve_fit


class BasicForestFire(object):
    def __init__(self, length, density=0):
        self.length = length

        self.forest = None
        self.reset(density=density)

    def reset(self, density):
        if density == 0:
            self.forest = np.zeros((self.length, self.length))
        else:
            total_spots = self.length ** 2
            nTrees = int(density * total_spots)
            arr = np.array([0] * (total_spots - nTrees) + [1] * nTrees)
            np.random.shuffle(arr)
            self.forest = arr.reshape((self.length, self.length))

    def start_fire(self, guarrantee_hit=False):
        if guarrantee_hit:
            xs, ys = np.nonzero(self.forest)
            index = np.random.randint(0, len(xs))
            x = xs[index]
            y = ys[index]
        else:
            x = np.random.randint(0, self.length)
            y = np.random.randint(0, self.length)
        if self.forest[x, y]:
            extended_forest = np.vstack((self.forest, self.forest[[0], :]))
            extended_forest = np.hstack((extended_forest,
                                         extended_forest[:, [0]]))
            labeled_forest, _ = label(extended_forest)
            fire_label = labeled_forest[x, y]
            labels = {fire_label}
            labels = labels.union(set(
                labeled_forest[0, labeled_forest[-1, :] == fire_label]))
            labels = labels.union(set(
                labeled_forest[-1, labeled_forest[0, :] == fire_label]))
            labels = labels.union(set(
                labeled_forest[labeled_forest[:, -1] == fire_label, 0]))
            labels = labels.union(set(
                labeled_forest[labeled_forest[:, 0] == fire_label, -1]))

            burned_matrix = np.in1d(labeled_forest[:-1, :-1], list(labels)
                                    ).reshape(self.forest.shape)
            self.forest[burned_matrix] = 0
            return np.nonzero(burned_matrix)


class ForestFire(BasicForestFire):
    def __init__(self, length, p, f, plot=False):
        self.plot_forest = plot
        self.p = p
        self.f = f
        super(ForestFire, self).__init__(length, density=0)

    def reset(self, density):
        super(ForestFire, self).reset(density=density)
        if self.plot_forest:
            self.setup_plot()

    def update(self, regrow=True):
        empty_sites = np.nonzero(self.forest == 0)
        burned_trees = None
        burned_ratio = 0
        regrown_burned_ratio = 0

        # Grow
        for x, y in zip(empty_sites[0], empty_sites[1]):
            if np.random.rand() < self.p:
                self.forest[x, y] = 1

        # Lightning
        if np.random.rand() < self.f:
            num_trees = np.sum(self.forest)
            density = num_trees / (self.length ** 2)
            burned_trees = self.start_fire()
            if burned_trees:
                burned_ratio = len(burned_trees[0]) / num_trees
            if burned_trees and regrow:
                regrown_forest = BasicForestFire(self.length, density=density)
                regrown_burned_trees = regrown_forest.start_fire(
                    guarrantee_hit=True)
                regrown_burned_ratio = len(regrown_burned_trees[0]) / num_trees

        if self.plot_forest:
            self.plot(burned_trees)

        return burned_ratio, regrown_burned_ratio

    def setup_plot(self):
        plt.ion()
        self.fig = plt.figure(figsize=[8, 8])

        subplot1 = self.fig.add_subplot(111)
        self.points = subplot1.plot([], [], 'g.', [], [], 'r.')
        subplot1.axis([1, self.length, 1, self.length])
        subplot1.set_title("p=%s, f=%s" % (self.p, self.f))

        # subplot2 = fig.add_subplot(122)
        # lines = subplot2.plot([], [], 'b-', [], [], 'r-', [], [], 'g-')
        # subplot2.set_ylim([0, 1])
        # subplot2.autoscale(True, 'x')
        # subplot2.legend([lines[0], lines[1], lines[2]], ['Susceptible',
        #                                                  'Infected',
        #                                                  'Recovered'])
        # subplot2.set_xlabel(r'Time')
        # subplot2.set_ylabel(r'Population density')
        # subplot2.set_title(r'd=%s, $\beta$=%s, $\gamma$=%s' %
        #                    (self.diffusion_rate, self.beta, self.gamma))

        self.fig.canvas.manager.show()

    def plot(self, burned_trees):
        inds = np.nonzero(self.forest)
        self.points[0].set_data(inds[0] + 1, inds[1] + 1)
        if burned_trees:
            burned_trees = [burned_trees[0] + 1, burned_trees[1] + 1]
        else:
            burned_trees = [tuple(), tuple()]
        self.points[1].set_data(burned_trees)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()


def plot_rank_freq(og_ratios, regrown_ratios, p, f, grid_size):
    og_ratios.sort(reverse=True)
    regrown_ratios.sort(reverse=True)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    og_rank = np.array(range(len(og_ratios))) / len(og_ratios)
    ax.plot(og_ratios, og_rank, 'bx', linewidth=1, markersize=1)
    regrown_rank = np.array(range(len(regrown_ratios))) / len(regrown_ratios)
    ax.plot(regrown_ratios, regrown_rank, 'rx', linewidth=1, markersize=1)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylabel('cCDF')
    ax.set_xlabel('Relative fire size')

    ax.set_title("grid-size=%s p=%s, f=%s" % (grid_size, p, f))

    ax.legend(("Original forest", "Regrown forest"))

    plt.show()


def display_fire_simulation(p, f):
    grid_size = 128
    forest_fire = ForestFire(grid_size, p, f, plot=True)
    for _ in range(2000):
        og_ratio, regrown_ratio = forest_fire.update()
        if og_ratio:
            plt.pause(0.1)
    print("Done simulating!")


def run_forest_fire_simulation(p,f):
    grid_size = 128
    forest_fire = ForestFire(grid_size, p, f)
    og_ratios = []
    regrown_ratios = []
    while len(regrown_ratios) < 1000:
        og_ratio, regrown_ratio = forest_fire.update()
        if og_ratio:
            og_ratios.append(og_ratio)
            regrown_ratios.append(regrown_ratio)
    print("Done simulating!")
    print(og_ratios)
    print(regrown_ratios)
    plot_rank_freq(og_ratios, regrown_ratios, p, f, grid_size)


def plot_and_fit_rank_freq(og_ratios, p, f, grid_size):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    og_ratios.sort(reverse=True)
    og_rank = (np.array(range(len(og_ratios))) + 1) / len(og_ratios)
    ax.plot(og_ratios, og_rank, 'bx', linewidth=1, markersize=1)

    # Calculate and plot linear fit
    cutoff = int(1 * len(og_ratios) / 3)
    coefficients = np.polyfit(np.log(og_ratios[cutoff:]),
                              np.log(og_rank[cutoff:]),
                              1)
    tao = 1 - coefficients[0]
    min_exponent = np.log10(min(og_ratios[cutoff:]))
    fit_x = np.logspace(min_exponent, 0)
    fit_y = np.exp(coefficients[0] * np.log(fit_x) + coefficients[1])
    ax.plot(fit_x, fit_y, 'r', linewidth=1, markersize=1)

    # Calculate and plot power law
    x_min = 1 / grid_size ** 2
    r = np.random.random(len(og_rank))
    pow_x = x_min * (1 - r) ** (-1 / 0.15)
    pow_x.sort()
    pow_x = pow_x[::-1]
    ax.plot(pow_x, og_rank, 'g', linewidth=1, markersize=1)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('cCDF')
    ax.set_xlabel('Relative fire size')
    ax.set_title("grid-size=%s p=%s, f=%s" % (grid_size, p, f))
    ax.legend(("Original forest", "Linear fit, tao=%.4f" % tao,
               "Power law, tao=1.15"))
    ax.set_xlim((min(og_ratios), 1))
    ax.set_ylim((min(og_rank), 1))

    plt.show()


def _calculate_og_ratios(p, f, grid_size, lightning_strikes):
    print("Simulating for grid size %s" % grid_size)
    forest_fire = ForestFire(grid_size, p, f)
    og_ratios = []
    while len(og_ratios) < lightning_strikes:
        # print(".", end='', flush=True)
        og_ratio, _ = forest_fire.update(regrow=False)
        if og_ratio:
            og_ratios.append(og_ratio)
            if len(og_ratios) % 100 == 0:
                print(len(og_ratios))
    print("Done simulating!")
    return og_ratios


def power_law_fit(p, f):
    grid_size = 128
    # print(og_ratios)
    og_ratios = _calculate_og_ratios(p, f, grid_size, lightning_strikes=1000)
    plot_and_fit_rank_freq(og_ratios, p, f, grid_size)


def _calculate_tao(og_ratios):
    og_ratios.sort(reverse=True)
    og_rank = (np.array(range(len(og_ratios))) + 1) / len(og_ratios)

    # Calculate and plot linear fit
    cutoff = int(1 * len(og_ratios) / 3)
    coefficients = np.polyfit(np.log(og_ratios[cutoff:]),
                              np.log(og_rank[cutoff:]),
                              1)
    tao = 1 - coefficients[0]
    return tao


def size_scaling(grid_sizes, p, f):
    p = 0.01
    f = 0.6
    taos = []
    for grid_size in grid_sizes:
        og_ratios = _calculate_og_ratios(p=p, f=f, grid_size=grid_size,
                                         lightning_strikes=1500)
        tao = _calculate_tao(og_ratios)
        taos.append(tao)
    print(taos)
    return taos


def plot_size_scaling(grid_sizes, taos, p, f):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(1 / grid_sizes, taos, 'bo')

    def fit_func(x, a, b, c):
        return a * x ** b + c

    params = curve_fit(fit_func, 1 / grid_sizes, taos)
    print(params)
    x = np.logspace(-4, np.log10(max(1 / grid_sizes)))
    y = fit_func(x, *params[0])
    ax.plot(x, y, 'r')
    lowest_tao = fit_func(0, *params[0])
    print("tau=%s for N->inf" % lowest_tao)
    ax.plot(x, x * 0 + lowest_tao, 'g--')

    ax.legend(("Simulation data", "Fitted curve",
               "tau=%.4f for N->inf" % lowest_tao))
    # ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('tao')
    ax.set_xlabel('1/N')
    ax.set_title("p=%s, f=%s" % (p, f))

    plt.show()


if __name__ == '__main__':
    p = 0.02
    f = 0.1

    # display_fire_simulation(p, f)

    run_forest_fire_simulation(p,f)
    # power_law_fit(p, f)

    # grid_sizes = np.array([8, 16, 32, 64, 128, 256, 512])
    # taos = size_scaling(grid_sizes, p, f)
    # taos = [1.5833218481405096, 1.3988072477793636, 1.2922491631353956,
    #         1.2289157115325862, 1.1980668038561126, 1.203210835311791,
    #         1.1857502968013702]
    # taos = [1.5551278856754545, 1.3814123163765233, 1.2890762094744472,
    #         1.2382936205456594, 1.2019862903648932, 1.2011314147454277,
    #         1.1699824583345639]
    # plot_size_scaling(grid_sizes, taos, p, f)
