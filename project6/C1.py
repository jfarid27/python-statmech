import math, random, pylab

class LevyFreeSimulation():

    def levy_free_path(self, xstart, xend, dtau, N):
        x = [xstart]
        for k in range(1, N):
            dtau_prime = (N - k) * dtau
            x_mean = (dtau_prime * x[k - 1] + dtau * xend) / \
                     (dtau + dtau_prime)
            sigma = math.sqrt(1.0 / (1.0 / dtau + 1.0 / dtau_prime))
            x.append(random.gauss(x_mean, sigma))
        return x

    def V(self, x, cubic, quartic):
        pot = x ** 2 / 2.0 + cubic * x ** 3 + quartic * x ** 4
        return pot

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
    
    def run(self, beta, cubic, quartic, N, n_steps=4000000, delta=1.0):
        dtau = float(beta) / N
        sigma = 1.0 / math.sqrt( 2.0 * math.tanh( beta / 2.0))
        x = [1.0] * N
        Ncut = N/2
        data = []
        for step in range(n_steps):
            x_new = self.levy_free_path(x[0], x[0], dtau, N)
            weight_x = math.exp(sum(-self.V(a, cubic, quartic) * dtau for a in x))
            weight_x_new = math.exp(sum(-self.V(a, cubic, quartic) * dtau for a in x_new))
            if (random.uniform(0.0, 1.0) < (weight_x_new/ weight_x) ):
                x = x_new[:]
                x = x[1:] + x[:1]
            if step % N == 0:
                k = random.randint(0, N - 1)
                data.append(x[k])
        return data, x

if __name__ == "__main__":
   beta = 20
   steps = 100
   cubic, quartic = (-0.0, 0.0)
   simulation = LevyFreeSimulation() 
   distribution, final_conf = simulation.run(beta, cubic, quartic, steps, n_steps=400000)
   simulation.plot_histogram(distribution, beta, steps)
