import math, random, pylab

def read_file(filename):
    list_x = []
    list_y = []
    with open(filename) as f:
        for line in f:
            x, y = line.split()
            list_x.append(float(x))
            list_y.append(float(y))
    f.close()
    return list_x, list_y

def rho_free(x, y, beta):    # free off-diagonal density matrix
    return math.exp(-(x - y) ** 2 / (2.0 * beta)) 

def V(x, cubic, quartic):
    return ((x ** 2) / 2.0) + (cubic * (x ** 3)) + (quartic * (x ** 4))

beta = 4.0
cubic = -1
quartic = 1
N = 8                                             # number of slices
dtau = beta / N
delta = 1.0                                       # maximum displacement on one slice
n_steps = 1000000                                 # number of Monte Carlo steps
x = [0.0] * N                                     # initial path
slices = []
for step in range(n_steps):
    k = random.randint(0, N - 1)                  # random slice
    knext, kprev = (k + 1) % N, (k - 1) % N       # next/previous slices
    x_new = x[k] + random.uniform(-delta, delta)  # new position at slice k
    old_weight  = (rho_free(x[knext], x[k], dtau) *
                   rho_free(x[k], x[kprev], dtau) *
                   math.exp(-0.5 * dtau * V(x[k], cubic, quartic) ))
    new_weight  = (rho_free(x[knext], x_new, dtau) *
                   rho_free(x_new, x[kprev], dtau) *
                   math.exp(-0.5 * dtau * V(x_new, cubic, quartic) ))
    if random.uniform(0.0, 1.0) < new_weight / old_weight:
        x[k] = x_new
    if (step % 10 == 0):
       slices.append(x[0]) 

x_dist_trotter, y_dist_trotter = read_file('data_anharm_matrixsquaring_beta4.0.dat')
pylab.xlim(-2.0, 2.0)
pylab.hist(slices, bins=100, normed = 'True')
pylab.plot(x_dist_trotter, y_dist_trotter, c='red', linewidth=2.0)
pylab.title('Anharmonic Trotter probability compared with Metropolis MCMC \
    \nhistogram for '+str(len(slices))+' samples', fontsize = 18)
pylab.xlabel('$x$', fontsize = 30)
pylab.ylabel('$\pi(x)$', fontsize = 30)
pylab.legend(['Trotter', 'MCMC'])
pylab.show()
