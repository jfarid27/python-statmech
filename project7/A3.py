import math, random, pylab

class TwoParticleBosonSimulation():

    def analyticDist(self, r, beta):
        """Returns theoretical probability for distinguishable
           particles.
        """ 
        sigma = math.sqrt(2.0) / math.sqrt(2.0 * math.tanh(beta / 2.0))
        prob = (math.sqrt(2.0 / math.pi) / sigma) * math.exp(- r ** 2 / 2.0 / sigma ** 2)
        return prob

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

    def runSeparationSimulation(self, beta=2.0, nsteps=1000000):
        """Returns an array of simulated separation distances for
           identical bosons in one dimensional harmonic trap
        """
        low_1, low_2 = self.levy_harmonic_path(2)
        x = {low_1:low_1, low_2:low_2}
        radiusData = []
        for step in xrange(nsteps):
            # move 1
            a = random.choice(x.keys())
            if a == x[a]:
                dummy = x.pop(a)
                a_new = self.levy_harmonic_path(1)[0]
                x[a_new] = a_new
            else:
                a_new, b_new = self.levy_harmonic_path(2)
                x = {a_new:b_new, b_new:a_new}
            # move 2
            (low1, high1), (low2, high2) = x.items()
            weight_old = self.rho_harm_1d(low1, high1, beta) * self.rho_harm_1d(low2, high2, beta)
            weight_new = self.rho_harm_1d(low1, high2, beta) * self.rho_harm_1d(low2, high1, beta)
            if random.uniform(0.0, 1.0) < weight_new / weight_old:
                x = {low1:high2, low2:high1}
            radiusData.append(abs(x.keys()[0] - x.keys()[1]))
        return radiusData

    def plotHistogramData(self, data, beta, N, show_theoretical=True, saveFig=False): 
        pylab.hist(data, normed=True, bins=150, label='Sampled')
        if (show_theoretical):
            list_x = [0.1 * a for a in range (0, 31)]
            list_y = [ self.analyticDist(x, beta) for x in list_x]
            pylab.plot(list_x, list_y, label='Analytic')
        pylab.legend()
        pylab.xlabel('$x$')
        pylab.ylabel('$\\pi(x)$ (normalized)')
        pylab.title('Two Particle Boson Simulation Pair Radius(beta=%s, N=%i)' % (beta, N))
        pylab.xlim(0, 2)
        if (saveFig):
            pylab.savefig('plot_A1_beta%s.png' % beta)
        pylab.show()

if (__name__ == '__main__'):
    beta = .1 
    nsteps = 1000000
    simulationInstance = TwoParticleBosonSimulation()
    data = simulationInstance.runSeparationSimulation(beta, nsteps)
    simulationInstance.plotHistogramData(data, beta, nsteps) 
