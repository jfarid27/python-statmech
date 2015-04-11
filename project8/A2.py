import random, math, pylab, os

class Ising2DSimulation():

    def generateConfiguration(self, L, T, readFromFile=False, filename="local_configuration.txt"):
        """Returns read configuration from filename or generates
           a new configuration.
        """
        if os.path.isfile(filename) and readFromFile:
            f = open(filename, 'r')
            S = []
            for line in f:
                S.append(int(line))
            f.close()
            print 'Starting from file', filename
        else:
            S = [random.choice([1, -1]) for k in range(L * L)]
            print 'Starting from a random configuration'
        return S

    def writeConfiguration(self, S, filename="local_configuration.txt"):
        """Writes given configuration to specified file
        """
        f = open(filename, 'w')
        for a in S:
           f.write(str(a) + '\n')
        f.close()

    def x_y(self, k, L):
        """Generates x and y coordinates from array index of 2D lattice
        """
        y = k // L
        x = k - y * L
        return x, y

    def plotConfiguration(self, S, L, saveFigure=False):
        """Plots configuration of given 2D Ising Lattice
        """
        conf = [[0 for x in range(L)] for y in range(L)]
        for k in range(L * L):
            x, y = self.x_y(k, L)
            conf[x][y] = S[k]
        pylab.imshow(conf, extent=[0, L, 0, L], interpolation='nearest')
        pylab.set_cmap('hot')
        pylab.title('Local_'+ str(T) + '_' + str(L))
        if (saveFigure):
            pylab.savefig('plot_A2_local_'+ str(T) + '_' + str(L)+ '.png')
        pylab.show()

    def energy(self, S, N, nbr):
        """Computes energy of given lattice with neighbor address
           dictionary
        """
        E = 0.0
        for k in range(N):
            E -=  S[k] * sum(S[nn] for nn in nbr[k])
        return 0.5 * E

    def simulateEnergy(self, S, L=6, T=2.0, nsteps=10000, \
        read=False, store=False, filename="local_configuration.txt", observeE=False):
        """Runs MCMC local simulation on 2D Ising lattice and returns energy
           observables
        """
        N = L * L
        nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                    (i // L) * L + (i - 1) % L, (i - L) % N) \
                                            for i in range(N)}
        beta = 1.0 / T
        Energy = self.energy(S, N, nbr)
        E = []
        for step in xrange(nsteps):
            k = random.randint(0, N - 1)
            delta_E = 2.0 * S[k] * sum(S[nn] for nn in nbr[k])
            if random.uniform(0.0, 1.0) < math.exp(-beta * delta_E):
                S[k] *= -1
                Energy += delta_E
            if(observeE):
                E.append(Energy)
        return E, S


if __name__ == "__main__":
    filename="A2_configuration.txt"
    T=1.0
    L=128
    nsteps=1000000
    instance = Ising2DSimulation()
    start_config = instance.generateConfiguration(L, T, True, filename)
    energy_samples, final_config = instance.simulateEnergy(start_config, L, T, nsteps)
    instance.plotConfiguration(final_config, L)
    instance.writeConfiguration(final_config, filename)
