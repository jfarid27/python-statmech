import random, math, os, pylab

class Wolff2DSimulation():

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

    def energy(self, S, N, nbr):
        """Computes energy of given lattice with neighbor address
           dictionary
        """
        E = 0.0
        for k in range(N):
            E -=  S[k] * sum(S[nn] for nn in nbr[k])
        return 0.5 * E

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

    def simulateEnergy(self, S, L=6, T=2.0, nsteps=10000, \
        observeE=False):
        """Runs MCMC local simulation on 2D Ising lattice and returns energy
           observables
        """

        N = L * L
        nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                    (i // L) * L + (i - 1) % L, (i - L) % N)
                                            for i in range(N)}
        p  = 1.0 - math.exp(-2.0 / T)
        mean_E = 0
        mean_E_per_N = 0
        mean_Esq = 0
        for step in range(nsteps):
            k = random.randint(0, N - 1)
            Pocket, Cluster = [k], [k]
            while Pocket != []:
                j = Pocket.pop()
                for l in nbr[j]:
                    if S[l] == S[j] and l not in Cluster \
                           and random.uniform(0.0, 1.0) < p:
                        Pocket.append(l)
                        Cluster.append(l)
            for j in Cluster:
                S[j] *= -1
            if (observeE):
                currEnergy = self.energy(S, N, nbr)
                mean_E_per_N += ( currEnergy/ ( float(nsteps) * float(N)))
                mean_E += ( currEnergy/ float(nsteps))
                mean_Esq += ( currEnergy ** 2 / float(nsteps) ) 
        return S, mean_E, mean_Esq, mean_E_per_N, \
                (mean_Esq - mean_E ** 2 ) / N / T ** 2 #cv calculation

if __name__ == "__main__":
    filename="B1_configuration.txt"
    T=2.27
    log2Ls = range(5)
    nsteps=100000
    instance = Wolff2DSimulation()

    def simmap(log2L):
        """The iterator to run simulations for a given L
        """
        L = 2 ** log2L
        start_config = instance.generateConfiguration(L, T)
        final_config, mean_E, mean_Esq, mean_E_per_N, cv \
            = instance.simulateEnergy(start_config, L, T, nsteps, True)
        print "Numbers for L = ", L
        print mean_E, mean_E_per_N, cv

    map(simmap, log2Ls)
