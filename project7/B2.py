import random, math, numpy, sys, os, pylab, mpl_toolkits.mplot3d

class MarkovChainHarmonicBosonSimulation():
    def levy_harmonic_path_3d(self, k, beta):
        x0 = tuple([random.gauss(0.0, 1.0 / math.sqrt(2.0 *
                    math.tanh(k * beta / 2.0))) for d in range(3)])
        x = [x0]
        for j in range(1, k):
            Upsilon_1 = 1.0 / math.tanh(beta) + 1.0 / \
                              math.tanh((k - j) * beta)
            Upsilon_2 = [x[j - 1][d] / math.sinh(beta) + x[0][d] /
                         math.sinh((k - j) * beta) for d in range(3)]
            x_mean = [Upsilon_2[d] / Upsilon_1 for d in range(3)]
            sigma = 1.0 / math.sqrt(Upsilon_1)
            dummy = [random.gauss(x_mean[d], sigma) for d in range(3)]
            x.append(tuple(dummy))
        return x

    def rho_harm_3d(self, x, xp, beta):
        Upsilon_1 = sum((x[d] + xp[d]) ** 2 / 4.0 *
                        math.tanh(beta / 2.0) for d in range(3))
        Upsilon_2 = sum((x[d] - xp[d]) ** 2 / 4.0 /
                        math.tanh(beta / 2.0) for d in range(3))
        return math.exp(- Upsilon_1 - Upsilon_2)

    def readConfiguration(self, beta, filename='data_boson_configuration.txt'):
        positions = {}
        if os.path.isfile(filename):
            f = open(filename, 'r')
            for line in f:
                a = line.split()
                positions[tuple([float(a[0]), float(a[1]), float(a[2])])] = \
                       tuple([float(a[3]), float(a[4]), float(a[5])])
            f.close()
            if len(positions) != N:
                sys.exit('ERROR in the input file.')
            print 'Starting from file', filename
        else:
            for k in range(N):
                a = self.levy_harmonic_path_3d(1, beta)
                positions[a[0]] = a[0]
            print 'Starting from a new configuration'
        return positions

    def writeConfiguration(self, positions, filename='data_boson_configuration.txt'):
        f = open(filename, 'w')
        for a in positions:
           b = positions[a]
           f.write(str(a[0]) + ' ' + str(a[1]) + ' ' + str(a[2]) + ' ' +
                   str(b[0]) + ' ' + str(b[1]) + ' ' + str(b[2]) + '\n')
        f.close() 

    def runSimulation(self, cycle_min=10, T_star=0.7, N=64, nsteps=10000, readConfiguration=False, writeConfiguration=False):
        beta = 1.0 / (T_star * N ** (1.0 / 3.0))
        # Initial condition
        if (readConfiguration):
            positions = self.readConfiguration(beta)
        else:
            positions = {}
            for k in range(N):
                a = self.levy_harmonic_path_3d(1, beta)
                positions[a[0]] = a[0]
        # Monte Carlo loop
        positionsX = []
        positionsXCycles = []
        for step in range(nsteps):
            # move 1: resample one permutation cycle
            boson_a = random.choice(positions.keys())
            perm_cycle = []
            while True:
                perm_cycle.append(boson_a)
                boson_b = positions.pop(boson_a)
                if boson_b == perm_cycle[0]:
                    break
                else:
                   boson_a = boson_b
            k = len(perm_cycle)
            perm_cycle = self.levy_harmonic_path_3d(k, beta)
            positions[perm_cycle[-1]] = perm_cycle[0]
            for k in range(len(perm_cycle) - 1):
                positions[perm_cycle[k]] = perm_cycle[k + 1]
            # move 2: exchange
            a_1 = random.choice(positions.keys())
            b_1 = positions.pop(a_1)
            a_2 = random.choice(positions.keys())
            b_2 = positions.pop(a_2)
            weight_new = self.rho_harm_3d(a_1, b_2, beta) * self.rho_harm_3d(a_2, b_1, beta)
            weight_old = self.rho_harm_3d(a_1, b_1, beta) * self.rho_harm_3d(a_2, b_2, beta)
            if random.uniform(0.0, 1.0) < weight_new / weight_old:
                positions[a_1] = b_2
                positions[a_2] = b_1
            else:
                positions[a_1] = b_1
                positions[a_2] = b_2
            positionsX.append(map(lambda boson: positions[boson][0], positions.keys())) 
            if (k > cycle_min):
                positionsXCycles.append(map(lambda boson: positions[boson][0], positions.keys())) 
        if (writeConfiguration):
            self.writeConfiguration(positions)
        return positions, positionsX, positionsXCycles
            

    def plot3dHistogramData(self, positions, T_star, N, saveFig=False): 
        fig = pylab.figure()
        ax = mpl_toolkits.mplot3d.axes3d.Axes3D(fig)
        ax.set_aspect('equal')
        n_colors = 10
        list_colors = pylab.cm.rainbow(numpy.linspace(0, 1, n_colors))[::-1]
        dict_colors = {}
        i_color = 0
        positions_copy = positions.copy()
        while positions_copy:
            x, y, z = [], [], []
            starting_boson = positions_copy.keys()[0]
            boson_old = starting_boson
            while True:
                x.append(boson_old[0])
                y.append(boson_old[1])
                z.append(boson_old[2])
                boson_new = positions_copy.pop(boson_old)
                if boson_new == starting_boson: break
                else: boson_old = boson_new
            len_cycle = len(x)
            if len_cycle > 2:
                x.append(x[0])
                y.append(y[0])
                z.append(z[0])
            if len_cycle in dict_colors:
                color = dict_colors[len_cycle]
                ax.plot(x, y, z, '+-', c=color, lw=0.75)
            else:
                color = list_colors[i_color]
                i_color = (i_color + 1) % n_colors
                dict_colors[len_cycle] = color
                ax.plot(x, y, z, '+-', c=color, label='k=%i' % len_cycle, lw=0.75)
        pylab.title(str(N) + ' bosons at T* = ' + str(T_star))
        pylab.legend()
        ax.set_xlabel('$x$', fontsize=16)
        ax.set_ylabel('$y$', fontsize=16)
        ax.set_zlabel('$z$', fontsize=16)
        xmax = 6.0
        ax.set_xlim3d([-xmax, xmax])
        ax.set_ylim3d([-xmax, xmax])
        ax.set_zlim3d([-xmax, xmax])
        if (saveFig):
            pylab.savefig('plot_boson_configuration.png')
        pylab.show()

    def plotHistogramData(self, data, T, N, saveFig=False): 
        beta = 1/float(T)
        pylab.hist(data, normed=True, bins=100, label='Sampled')
        pylab.legend()
        pylab.xlabel('$x$')
        pylab.ylabel('$\\pi(x)$ (normalized)')
        pylab.title('Boson Simulation X Positions(beta=%s, N=%i)' % (beta, N))
        pylab.xlim(-3, 3)
        if (saveFig):
            pylab.savefig('plot_A1_beta%s.png' % beta)
        pylab.show()

if (__name__ == '__main__'):
    T = 0.8 
    nsteps = 10000
    nparticles = 64
    simulationInstance = MarkovChainHarmonicBosonSimulation()
    data, xPositions, xCycles = simulationInstance.runSimulation(T, nparticles, nsteps)
    simulationInstance.plotHistogramData(xPositions, T, nparticles) 
    simulationInstance.plotHistogramData(xCycles, T, nparticles) 
