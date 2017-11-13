import matplotlib.pyplot as plt
import numpy as np


class Lattice(object):
    def __init__(self, length, num_total_individuals, num_infected, beta,
                 gamma, diffusion_rate):
        self.length = length
        self.beta = beta
        self.gamma = gamma
        self.diffusion_rate = diffusion_rate

        self.num_individuals = num_total_individuals
        self.num_initially_infected = num_infected

        self.num_infected = 0
        self.num_susceptible = 0
        self.num_recovered = 0

        self.lattice = self._get_clear_lattice()

        self.reset()

    def reset(self):
        self.num_infected = self.num_initially_infected
        self.num_susceptible = self.num_individuals - self.num_infected
        self.num_recovered = 0

        self.lattice = self._get_clear_lattice()

        for i in range(self.num_individuals):
            x = np.random.randint(0, self.length)
            y = np.random.randint(0, self.length)
            if i < self.num_infected:
                self.lattice[x, y, 1] = 1
            else:
                self.lattice[x, y, 0] += 1

    def run_simulation(self, optimize=True, plot_during=False,
                       plot_final=False):
        i = 0
        susceptibles = []
        infected = []
        recovered = []

        if plot_during:
            fig1, ax1 = plt.subplots()

        while self.num_infected:
            self.random_walk(optimize=optimize)
            self.check_infection()
            if plot_during and i % 20 == 0:
                self.plot(ax1)
            if i % 10000 == 0 and i != 0:
                print('        simulation has run %s time steps' % i)
                print('        infected: %s, recovered: %s' %
                      (self.num_infected, self.num_recovered))
            susceptibles.append(self.num_susceptible)
            infected.append(self.num_infected)
            recovered.append(self.num_recovered)
            i += 1

        if plot_final:
            fig, ax = plt.subplots()
            ax.plot(susceptibles, 'b')
            ax.plot(recovered, 'g')
            ax.plot(infected, 'r')
            plt.show()
        return self.num_recovered / self.num_individuals

    def plot(self, ax):
        ax.clear()
        ax.axis([0, self.length - 1, 0, self.length - 1])
        susceptible_inds = np.nonzero(self.lattice[:, :, 0])
        ax.scatter(susceptible_inds[0], susceptible_inds[1], c='b')
        infected_inds = np.nonzero(self.lattice[:, :, 1])
        ax.scatter(infected_inds[0], infected_inds[1], c='r')
        recovered_inds = np.nonzero(self.lattice[:, :, 2])
        ax.scatter(recovered_inds[0], recovered_inds[1], c='g')
        plt.pause(0.01)

    def random_walk(self, optimize):
        new_lattice = self._get_clear_lattice()

        if optimize:
            new_lattice[:, :, 2] = self.lattice[:, :, 2]
            nonzero_indices = np.nonzero(np.sum(self.lattice[:, :, :2], 2))
        else:
            nonzero_indices = np.nonzero(np.sum(self.lattice, 2))

        for x, y in zip(nonzero_indices[0], nonzero_indices[1]):

            for _ in range(self.lattice[x, y, 0]):
                if np.random.rand() < self.diffusion_rate:
                    new_x, new_y = self._get_next_position(x, y)
                else:
                    new_x, new_y = x, y
                new_lattice[new_x, new_y, 0] += 1

            for _ in range(self.lattice[x, y, 1]):
                if np.random.rand() < self.diffusion_rate:
                    new_x, new_y = self._get_next_position(x, y)
                else:
                    new_x, new_y = x, y
                new_lattice[new_x, new_y, 1] += 1

            if not optimize:
                for _ in range(self.lattice[x, y, 2]):
                    if np.random.rand() < self.diffusion_rate:
                        new_x, new_y = self._get_next_position(x, y)
                    else:
                        new_x, new_y = x, y
                    new_lattice[new_x, new_y, 2] += 1

        self.lattice = new_lattice

    def check_infection(self):
        nonzero_indices = np.nonzero(self.lattice[:, :, 1])
        if len(nonzero_indices[0]) == 0:
            assert self.num_infected == 0
            return
        for x, y in zip(nonzero_indices[0], nonzero_indices[1]):
            num_infected = self.lattice[x, y, 1]
            for i in range(num_infected):
                if self.lattice[x, y, 0] and np.random.rand() < self.beta:
                    # Infect all susceptible individuals at the site
                    self.lattice[x, y, 1] += self.lattice[x, y, 0]
                    self.num_infected += self.lattice[x, y, 0]
                    self.num_susceptible -= self.lattice[x, y, 0]
                    self.lattice[x, y, 0] = 0
                if np.random.rand() < self.gamma:
                    # Recover (remove from lattice)
                    self.lattice[x, y, 1] -= 1
                    self.lattice[x, y, 2] += 1
                    self.num_infected -= 1
                    self.num_recovered += 1

    def _get_next_position(self, x, y):
        r = np.random.rand()
        if r < 1 / 4:
            x -= 1
        if r < 0.5:
            x += 1
        if r < 0.75:
            y -= 1
        else:
            y += 1

        # We don't need to check the lower bound since -1 is an allowed index
        x = 0 if x == self.length else x
        y = 0 if y == self.length else y

        return x, y

    def _get_clear_lattice(self):
        return np.zeros((self.length, self.length, 3), dtype=np.int8)
