import lattice


l = lattice.Lattice(length=100,
                    num_total_individuals=1000,
                    num_infected=10,
                    beta=0.5,
                    diffusion_rate=0.3,
                    gamma=0.005)
l.initialize_individuals()
l.run_simulation(plot_during=False, plot_final=False)
