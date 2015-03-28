import random, math, pylab

class QuantumSimulation:
    """
    Quantum simulation class for a harmonic oscillator at energy level n.
    """

    def __init__(self, beta):
        """
        Initialize the simulation with given inverse temperature constant
        beta
        """
        self.beta = beta
    
    def psi_n_square(self, n, x):
        """
        Recursive wavefunction calculator.
        """
        if n == -1:
            return 0.0
        else:
            psi = [math.exp(-x ** 2 / 2.0) / math.pi ** 0.25]
            psi.append(math.sqrt(2.0) * x * psi[0])
            for k in range(2, n + 1):
                psi.append(math.sqrt(2.0 / k) * x * psi[k - 1] -
                           math.sqrt((k - 1.0) / k) * psi[k - 2])
            return psi[n] ** 2

    def wavefunction(self, *args):
        """ 
        Harmonic oscillator wavefunction at energy level n.
        Here the first given arg is energy level and second is position.
        """
        E_n = args[0] + .5
        return self.psi_n_square(*args) * math.exp(-self.beta*E_n)

    def x_density(self, *params):
        """
        Computes probability of being at position x.
        """
        return self.wavefunction(*params) * self.wavefunction(*params)

    def analytical_density(self, x):
        """
        Analytical density function of quantum simulator
        """
        return math.sqrt(math.tanh(float(self.beta)/2) / math.pi) \
           * math.exp( -( x ** 2) * math.tanh(float(self.beta)/2) )

    def classical_density(self, x):
        """
        Classical density function of quantum simulator
        """
        return math.sqrt(float(self.beta)/ (2*math.pi)) \
            * math.exp(-self.beta * (x ** 2) / float(2) )

    def markov_metropolis(self, x=0.0, n=0, delta=.05, \
        transitions=100000, hist=False):
        """Returns stop position, stop energy level, and optional lists of
        step positions and step energy levels.

        A function to compute the position of a particle using markov chain
        metropolis algorithm. Even transition numbers update energy level
        while odd transition numbers update positions.
        """
        position_values = []
        energy_values = []
        for k in range(transitions):
            if (k % 2 == 0):
              #Update energy levels
                energy_options = [n-1, n+1]
                transition_n = random.choice(energy_options)
                if random.uniform(0.0, 1.0) <  \
                self.x_density(transition_n, x) / self.x_density(n, x): 
                    n = transition_n 
            else:
              #Update position
                x_new = x + random.uniform(-delta, delta)
                if random.uniform(0.0, 1.0) <  \
                self.x_density(n, x_new) / self.x_density(n, x): 
                    x = x_new 
            if hist: 
                position_values.append(x) 
                energy_values.append(n)
        return x,n, position_values, energy_values

    def output_sim_plot(self):
        """
        Function to output a markov simulation plot.
        """
        end_x, end_n, positions, energies =  \
            self.markov_metropolis(hist=True, transitions=1000000)
        pylab.hist(positions, 100, normed = 'True')
        x = [a / 10.0 for a in range(-50, 51)]
        y = [self.analytical_density(a) for a in x]
        z = [self.classical_density(a) for a in x]
        pylab.plot(x, y, c='red', linewidth=2.0)
        pylab.plot(x, z, c='green', linewidth=2.0)
        pylab.title('Theoretical Gaussian distribution $\pi(x)$ and \
            \nhistogram for '+str(len(positions))+' samples, beta='+str(self.beta), fontsize = 18)
        pylab.xlabel('$x$', fontsize = 30)
        pylab.ylabel('$\pi(x)$', fontsize = 30)
        pylab.savefig('markov_sim_beta_'+str(self.beta)+'.png')


if __name__ == "__main__":
    beta = 5
    QuantumSimulation(beta=beta).output_sim_plot()
