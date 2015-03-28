import math, random, pylab

class LevyMCMCSimulation():
    def plot_configuration(self, configuration, dtau, beta, N):
        ys = [n * dtau for n in range(0,N)]
        xs = [x for x in configuration]
        pylab.plot(xs, ys, label='Path')
        pylab.legend()
        pylab.xlabel('$x$')
        pylab.ylabel('$Beta$')
        pylab.title('naive_harmonic_path (beta=%s, N=%i)' % (beta, N))
        pylab.xlim(-2, 2)
        pylab.show()
    
    def levy_harmonic_path(self,xstart, xend, dtau, N):
        x = [xstart]
        for k in range(1, N):
            dtau_prime = (N - k) * dtau
            Ups1 = 1.0 / math.tanh(dtau) + \
                   1.0 / math.tanh(dtau_prime)
            Ups2 = x[k - 1] / math.sinh(dtau) + \
                   xend / math.sinh(dtau_prime)
            x.append(random.gauss(Ups2 / Ups1, \
                   1.0 / math.sqrt(Ups1)))
        return x
    
    def plot_histogram(self, data, beta, N):
        pylab.hist(data, normed=True, bins=100, label='QMC')
        list_x = [0.1 * a for a in range (-30, 31)]
        list_y = [math.sqrt(math.tanh(beta / 2.0)) / math.sqrt(math.pi) * \
                  math.exp(-x ** 2 * math.tanh(beta / 2.0)) for x in list_x]
        pylab.plot(list_x, list_y, label='analytic')
        pylab.legend()
        pylab.xlabel('$x$')
        pylab.ylabel('$\\pi(x)$ (normalized)')
        pylab.title('naive_harmonic_path (beta=%s, N=%i)' % (beta, N))
        pylab.xlim(-2, 2)
        pylab.show()

    def rho_free(self, x, y, beta):
        return math.exp(-(x - y) ** 2 / (2.0 * beta))
    
    def run(self, beta, N, n_steps=4000000, delta=1.0):
        dtau = float(beta) / N
        x = [5.0] * N
        Ncut = N/2
        data = []
        for step in range(n_steps):
            x = self.levy_harmonic_path(x[0], x[0], dtau, N)
            x = x[Ncut:] + x[:Ncut]
            if step % N == 0:
                k = random.randint(0, N - 1)
                data.append(x[k])
        return data, x

if __name__ == "__main__":
   beta = 20
   steps = 80
   simulation = LevyMCMCSimulation() 
   distribution, final_conf = simulation.run(beta, steps, n_steps=400000)
   simulation.plot_histogram(distribution, beta, steps)
