import random

import numpy as np
import matplotlib.pyplot as plt


class Lattice(object):
    def __init__(self, length, num_total_individuals, num_infected, beta,
                 gamma, diffusion_rate):
        self.length = length
        self.beta = beta
        self.gamma = gamma
        self.diffusion_rate = diffusion_rate

        self.num_individuals = num_total_individuals
        self.num_infected = num_infected
        self.num_susceptible = self.num_individuals - self.num_infected
        self.num_recovered = 0

        self.lattice = self._get_clear_lattice()
        self.fig, self.ax = plt.subplots()

    def initialize_individuals(self):
        for i in range(self.num_individuals):
            x = random.randint(0, self.length-1)
            y = random.randint(0, self.length-1)
            if i < self.num_infected:
                self.lattice[x, y, 1] = 1
            else:
                self.lattice[x,y,0] += 1

    def run_simulation(self, plot_during, plot_final=True):
        i = 0
        susceptibles = []
        infected = []
        recovered = []
        while self.num_infected:

            self.random_walk()
            self.check_infection()
            if plot_during and i % 10 == 0:
                self.plot()
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
        else:
            print(self.num_recovered/self.num_individuals)

    def plot(self):
        self.ax.clear()
        self.ax.axis([0, self.length-1, 0, self.length-1])
        susceptible_inds = np.nonzero(self.lattice[:,:,0])
        self.ax.scatter(susceptible_inds[0], susceptible_inds[1], c='b')
        infected_inds = np.nonzero(self.lattice[:,:,1])
        self.ax.scatter(infected_inds[0], infected_inds[1], c='r')
        recovered_inds = np.nonzero(self.lattice[:,:,2])
        self.ax.scatter(recovered_inds[0], recovered_inds[1], c='g')
        plt.pause(0.01)

    def random_walk(self):
        new_lattice = self._get_clear_lattice()
        for x in range(self.length):
            for y in range(self.length):
                for _ in range(self.lattice[x,y,0]):
                    new_x, new_y = self._get_next_position(x,y)
                    new_lattice[new_x,new_y,0] += 1
                for _ in range(self.lattice[x,y,1]):
                    new_x, new_y = self._get_next_position(x,y)
                    new_lattice[new_x,new_y,1] += 1
                for _ in range(self.lattice[x, y, 2]):
                    new_x, new_y = self._get_next_position(x, y)
                    new_lattice[new_x, new_y, 2] += 1
        self.lattice = new_lattice

    def check_infection(self):
        for x in range(self.length):
            for y in range(self.length):
                num_infected = self.lattice[x,y,1]
                for i in range(num_infected):
                    if random.random() < self.beta:
                        # Infect all susceptible individuals at the site
                        self.lattice[x,y,1] += self.lattice[x,y,0]
                        self.num_infected += self.lattice[x,y,0]
                        self.num_susceptible -= self.lattice[x,y,0]
                        self.lattice[x,y,0] = 0
                    if random.random() < self.gamma:
                        # Recover (remove from lattice)
                        self.lattice[x,y,1] -= 1
                        self.lattice[x,y,2] += 1
                        self.num_infected -= 1
                        self.num_recovered += 1

    def _get_next_position(self, x, y):
        if random.random() > self.diffusion_rate:
            return x, y
        r = random.random()
        if r < 1/4:
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
        return np.zeros((self.length,self.length,3), dtype=np.int8)
