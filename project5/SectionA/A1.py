import random, math

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

    def markov_metropolis(self, probability, delta=.05, start=0.0, steps=100000):
        """
        A function to compute the position of a particle using markov chain
        metropolis algorithm.
        """
        for k in range(steps):
            x_new = x + random.uniform(-delta, delta)
            if random.uniform(0.0, 1.0) <  \
              prob_dist(x_new) / prob_dist(x): 
                x = x_new 
        return x

if __name__ == "__main__":
    print "HELLO"
