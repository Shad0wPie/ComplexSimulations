import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.measurements import label


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

    def update(self):
        empty_sites = np.nonzero(self.forest == 0)
        burned_trees = None
        regrown_burned_trees = None
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
                regrown_forest = BasicForestFire(self.length, density=density)
                regrown_burned_trees = regrown_forest.start_fire(
                    guarrantee_hit=True)
                burned_ratio = len(burned_trees[0]) / num_trees
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


def run_forest_fire_simulation():
    p = 0.01
    f = 0.6
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


if __name__ == '__main__':
    run_forest_fire_simulation()
