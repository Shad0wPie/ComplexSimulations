import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.measurements import label


class ForestFire(object):
    def __init__(self, length, p, f):
        self.length = length
        self.p = p
        self.f = f

        self.forest = None
        self.reset()

    def reset(self):
        self.forest = np.matlib.zeros((self.length, self.length))
        self.setup_plot()

    def update(self):
        empty_sites = np.nonzero(self.forest == 0)
        burned_trees = None
        for x, y in zip(empty_sites[0], empty_sites[1]):
            if np.random.rand() < self.p:
                self.forest[x, y] = 1
        if np.random.rand() < self.f:
            burned_trees = self.start_fire()
        self.plot(burned_trees)

    def start_fire(self):
        x = np.random.randint(0, self.length)
        y = np.random.randint(0, self.length)
        if self.forest[x, y]:
            extended_forest = np.vstack((self.forest, self.forest[0, :]))
            extended_forest = np.hstack((extended_forest,
                                         extended_forest[:, 0]))
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


def run_forest_fire_simulation():
    forest_fire = ForestFire(128, 0.01, 0.6)
    while True:
        forest_fire.update()


if __name__ == '__main__':
    run_forest_fire_simulation()
