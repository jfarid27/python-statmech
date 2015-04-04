import math, random, pylab

class TwoParticleDistinguishableSimulation():

    def analyticDist(self, x, beta):
        sigma = 1.0 / math.sqrt(2.0 * math.tanh(beta / 2.0))
        return math.exp(-x ** 2 / (2.0 * sigma ** 2)) / math.sqrt(2.0 * math.pi) / sigma

    def levy_harmonic_path(self, k):
        x = [random.gauss(0.0, 1.0 / math.sqrt(2.0 * math.tanh(k * beta / 2.0)))]
        if k == 2:
            Ups1 = 2.0 / math.tanh(beta)
            Ups2 = 2.0 * x[0] / math.sinh(beta)
            x.append(random.gauss(Ups2 / Ups1, 1.0 / math.sqrt(Ups1)))
        return x[:]

    def runSimulation(self, beta=2.0, nsteps=1000000):
        low = self.levy_harmonic_path(2)
        high = low[:]
        data = []
        for step in xrange(nsteps):
            k = random.choice([0, 1])
            low[k] = self.levy_harmonic_path(1)[0]
            high[k] = low[k]
            data.append(low[k])
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
        pylab.title('Two Particle Distinguishable Simulation(beta=%s, N=%i)' % (beta, N))
        pylab.xlim(-2, 2)
        if (saveFig):
            pylab.savefig('plot_A1_beta%s.png' % beta)
        pylab.show()

if (__name__ == '__main__'):
    beta = 2.0 
    nsteps = 100000
    simulationInstance = TwoParticleDistinguishableSimulation()
    data = simulationInstance.runSimulation(beta, nsteps)
    simulationInstance.plotHistogramData(data, beta, nsteps) 
