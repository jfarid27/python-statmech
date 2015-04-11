import random, math, pylab

class Heatbath2DIsing():

    def generateFixedConfiguration(self, sign, L):
        """Returns fixed configuration of size L^2
        """
        return [sign] * L * L

    def show_spins(self, S0, S1, L, label, save=False):
        N = L * L
        pylab.set_cmap('hot')
        conf0 = [[0 for x in range(L)] for y in range(L)]
        conf1 = [[0 for x in range(L)] for y in range(L)]
        for k in range(N):
            y = k // L
            x = k - y * L
            conf0[x][y] = S0[k]
            conf1[x][y] = S1[k]
        pylab.subplot(1, 2, 1)
        pylab.imshow(conf0, extent=[0, L, 0, L], \
            interpolation='nearest')
        pylab.title('S0 ' + label)
        pylab.subplot(1, 2, 2)
        pylab.imshow(conf1, extent=[0, L, 0, L], \
            interpolation='nearest')
        pylab.title('S1 ' + label)
        pylab.tight_layout()
        if (save):
            pylab.savefig('plot_' + label + '.png')
        pylab.close()

    def simulatePair(self, S0, S1, L=6, T=2.0):
        N = L * L
        step = 0
        nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                    (i // L) * L + (i - 1) % L, (i - L) % N) \
                                            for i in range(N)}
        beta = 1.0 / T
        while True:
            step += 1
            k = random.randint(0, N - 1)
            Upsilon = random.uniform(0.0, 1.0)
            h = sum(S0[nn] for nn in nbr[k])
            S0[k] = -1
            if Upsilon < 1.0 / (1.0 + math.exp(-2.0 * beta * h)):
                S0[k] = 1
            h = sum(S1[nn] for nn in nbr[k])
            S1[k] = -1
            if Upsilon < 1.0 / (1.0 + math.exp(-2.0 * beta * h)):
                S1[k] = 1
            if step % N == 0:
                n_diff = sum(abs(S0[i] - S1[i]) for i in range(N))
                if n_diff == 0:
                    t_coup = step / N
                    break
        return S0, S1, t_coup

    def multiSimPair(self, simulations, L=6, T=3.0):
        couplingTimes = []
        for k in range(simulations):
            S0 = self.generateFixedConfiguration(-1, L)
            S1 = self.generateFixedConfiguration(1, L)
            finalS1, finalS2, couplingTime = self.simulatePair(S0,S1,L,T)
            couplingTimes.append(couplingTime)
        print "Average coupling time for T=%f:" % T
        print sum(couplingTimes) / float(simulations)

if __name__ == "__main__":
    L=32
    T=[5.0, 4.0, 3.0, 2.5, 2.4]
    instance = Heatbath2DIsing()
    for t in T:
        instance.multiSimPair(20, L, t)
