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

        middle = int(self.length / 2)
        for i in range(self.num_individuals):
            if i < self.num_infected:
                self.lattice[middle, middle, 1] += 1
            else:
                x = np.random.randint(0, self.length)
                y = np.random.randint(0, self.length)
                self.lattice[x, y, 0] += 1

    def run_simulation(self, optimize=True, plot=False, plotstep=20):
        i = 0
        susceptibles = []
        infected = []
        recovered = []

        if plot:
            plt.ion()
            fig = plt.figure(figsize=[8, 4])

            subplot1 = fig.add_subplot(121)
            points = subplot1.plot([], [], 'b.', [], [], 'r.', [], [], 'g.')
            subplot1.axis([1, self.length, 1, self.length])

            subplot2 = fig.add_subplot(122)
            lines = subplot2.plot([], [], 'b-', [], [], 'r-', [], [], 'g-')
            subplot2.set_ylim([0, 1])
            subplot2.autoscale(True, 'x')
            subplot2.legend([lines[0], lines[1], lines[2]], ['Susceptible',
                                                             'Infected',
                                                             'Recovered'])
            subplot2.set_xlabel(r'Time')
            subplot2.set_ylabel(r'Population density')
            subplot2.set_title(r'd=%s, $\beta$=%s, $\gamma$=%s' %
                               (self.diffusion_rate, self.beta, self.gamma))

            fig.canvas.manager.show()

        while self.num_infected:
            if self.num_susceptible or (not optimize):
                self.random_walk(optimize=optimize)
                self.check_infection()
                if plot and i % plotstep == 0:
                    self.plot(points, lines, i)
                    subplot2.relim()
                    subplot2.autoscale_view(True, True, True)
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                if i % 10000 == 0 and i != 0:
                    print('        simulation has run %s time steps' % i)
                    print('        infected: %s, recovered: %s' %
                          (self.num_infected, self.num_recovered))
            else:
                self.num_recovered += self.num_infected
                self.num_infected = 0
            susceptibles.append(self.num_susceptible)
            infected.append(self.num_infected)
            recovered.append(self.num_recovered)
            i += 1
        print('        simulation completed after %s time steps' % i)
        if plot:
            plt.show(block=True)
        return self.num_recovered / self.num_individuals

    def plot(self, points, lines, iteration_number):
        for i in range(3):
            inds = np.nonzero(self.lattice[:, :, i])
            points[i].set_data(inds[0]+1, inds[1]+1)
            count = self.num_susceptible if i == 0 else self.num_infected if \
                i == 1 else self.num_recovered
            count /= self.num_individuals
            xdata = lines[i].get_xdata()
            lines[i].set_xdata(np.append(xdata, iteration_number))
            lines[i].set_ydata(np.append(lines[i].get_ydata(), count))

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
            assert self.num_infected == 0, '%s, %s, %s' % (self.num_infected,
                                                           self.num_recovered,
                                                           self.num_susceptible)
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
        elif r < 0.5:
            x += 1
        elif r < 0.75:
            y -= 1
        else:
            y += 1

        # We don't need to check the lower bound since -1 is an allowed index
        x = 0 if x == self.length else x
        y = 0 if y == self.length else y

        return x, y

    def _get_clear_lattice(self):
        return np.zeros((self.length, self.length, 3), dtype=np.int8)
