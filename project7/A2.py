import math, random, pylab

class TwoParticleBosonSimulation():

    def z(self, beta):
        return 1.0 / (1.0 - math.exp(- beta))

    def analyticDist(self, x, beta):
        pi_x_1 = math.sqrt(math.tanh(beta / 2.0)) / math.sqrt(math.pi) *\
                 math.exp(-x ** 2 * math.tanh(beta / 2.0))
        pi_x_2 = math.sqrt(math.tanh(beta)) / math.sqrt(math.pi) *\
                 math.exp(-x ** 2 * math.tanh(beta))
        weight_1 = self.z(beta) ** 2 / (self.z(beta) ** 2 + self.z(2.0 * beta))
        weight_2 = self.z(2.0 * beta) / (self.z(beta) ** 2 + self.z(2.0 * beta))
        pi_x = pi_x_1 * weight_1 + pi_x_2 * weight_2
        return pi_x

    def levy_harmonic_path(self, k):
        x = [random.gauss(0.0, 1.0 / math.sqrt(2.0 * math.tanh(k * beta / 2.0)))]
        if k == 2:
            Ups1 = 2.0 / math.tanh(beta)
            Ups2 = 2.0 * x[0] / math.sinh(beta)
            x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
        return x[:]

    def rho_harm_1d(self, x, xp, beta):
        Upsilon_1 = (x + xp) ** 2 / 4.0 * math.tanh(beta / 2.0)
        Upsilon_2 = (x - xp) ** 2 / 4.0 / math.tanh(beta / 2.0)
        return math.exp(- Upsilon_1 - Upsilon_2)

    def runSimulation(self, beta=2.0, nsteps=1000000):
        low = self.levy_harmonic_path(2)
        high = low[:]
        data = []
        for step in xrange(nsteps):
            # move 1
            if low[0] == high[0]:
                k = random.choice([0, 1])
                low[k] = self.levy_harmonic_path(1)[0]
                high[k] = low[k]
            else:
                low[0], low[1] = self.levy_harmonic_path(2)
                high[1] = low[0]
                high[0] = low[1]
            data += low[:]
            # move 2
            weight_old = (self.rho_harm_1d(low[0], high[0], beta) *
                          self.rho_harm_1d(low[1], high[1], beta))
            weight_new = (self.rho_harm_1d(low[0], high[1], beta) *
                          self.rho_harm_1d(low[1], high[0], beta))
            if random.uniform(0.0, 1.0) < weight_new / weight_old:
                high[0], high[1] = high[1], high[0]
        return data

    def plotHistogramData(self, data, beta, N, show_theoretical=True, saveFig=False): 
        pylab.hist(data, normed=True, bins=100, label='Sampled')
        if (show_theoretical):
            list_x = [0.1 * a for a in range (-30, 31)]
            list_y = [ self.analyticDist(x, beta) for x in list_x]
            pylab.plot(list_x, list_y, label='Analytic')
        pylab.legend()
        pylab.xlabel('$x$')
        pylab.ylabel('$\\pi(x)$ (normalized)')
        pylab.title('Two Particle Boson Simulation(beta=%s, N=%i)' % (beta, N))
        pylab.xlim(-2, 2)
        if (saveFig):
            pylab.savefig('plot_A1_beta%s.png' % beta)
        pylab.show()

if (__name__ == '__main__'):
    beta = 2.0 
    nsteps = 100000
    simulationInstance = TwoParticleBosonSimulation()
    data = simulationInstance.runSimulation(beta, nsteps)
    simulationInstance.plotHistogramData(data, beta, nsteps) 
