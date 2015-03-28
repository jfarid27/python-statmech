import random, math, pylab

class GroundStateSimulation():
    """
    Simulation class for the ground state wavefunction in harmonic
    oscillator potential.
    """
    
    def wavefunction(self, x):
        """
        Ground state wavefunction for a particle in a harmonic oscillator
        potential.
        """
        return ( float(1) / math.pi**(float(1)/4)) * math.exp( x**2 / float(-2))
    def x_density_function(self, x):
        """
        Computes probability of being at position x.
        """
        return self.wavefunction(x) * self.wavefunction(x)

    def markov_metropolis(self, x=0, delta=.05, start=0.0, \
        hist=False, steps=100000):
        """
        A function to compute the position of a particle using markov chain
        metropolis algorithm.
        """
        position_values = []
        for k in range(steps):
            x_new = x + random.uniform(-delta, delta)
            if random.uniform(0.0, 1.0) <  \
            self.x_density_function(x_new) / self.x_density_function(x): 
                x = x_new 
            if hist: position_values.append(x) 
        return x, position_values

if __name__ == "__main__":
    simulation_instance = GroundStateSimulation()
    endpoint, data =  simulation_instance.markov_metropolis(hist=True) 
    pylab.hist(data, 100, normed = 'True')
    x = [a / 10.0 for a in range(-50, 51)]
    y = [simulation_instance.x_density_function(a) for a in x]
    pylab.plot(x, y, c='red', linewidth=2.0)
    pylab.title('Theoretical Gaussian distribution $\pi(x)$ and \
        \nhistogram for '+str(len(data))+' samples', fontsize = 18)
    pylab.xlabel('$x$', fontsize = 30)
    pylab.ylabel('$\pi(x)$', fontsize = 30)
    pylab.show()
