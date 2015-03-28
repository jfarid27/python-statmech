import math, numpy, pylab

def V(x, cubic, quartic):
    return ((x ** 2) / 2.0) + (cubic * (x ** 3)) + (quartic * (x ** 4))

# Free off-diagonal density matrix
def rho_free(x, xp, beta):
    return (math.exp(-(x - xp) ** 2 / (2.0 * beta)) /
            math.sqrt(2.0 * math.pi * beta))

# Harmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_harmonic_trotter(grid, beta):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp(-0.5 * beta * 0.5 * (x ** 2 + xp ** 2)) \
                         for x in grid] for xp in grid])

# Anharmonic density matrix in the Trotter approximation (returns the full matrix)
def rho_anharmonic_trotter(grid, beta, cubic, quartic):
    return numpy.array([[rho_free(x, xp, beta) * \
                         numpy.exp(-0.5 * beta * (V(x, cubic, quartic) + V(xp, cubic, quartic))) \
                         for x in grid] for xp in grid])

x_max = 5.0                              # the x range is [-x_max,+x_max]
nx = 100
dx = 1.0 * x_max / (nx - 1)
x = [i * dx for i in range(-(nx - 1) / 2, nx / 2 + 1)]
beta_tmp = 2.0 ** (-5)                   # initial value of beta (power of 2)
beta_tmp_an = 2.0 ** (-5)                   # initial value of beta (power of 2)
beta     = 2.0 ** 2                      # actual value of beta (power of 2)
cubic = -1
quartic = 1

rho = rho_harmonic_trotter(x, beta_tmp)  # density matrix at initial beta
while beta_tmp < beta:
    rho = numpy.dot(rho, rho)
    rho *= dx
    beta_tmp *= 2.0

rho_an = rho_anharmonic_trotter(x, beta_tmp_an, cubic, quartic)  # density matrix at initial beta
while beta_tmp_an < beta:
    rho_an = numpy.dot(rho_an, rho_an)
    rho_an *= dx
    beta_tmp_an *= 2.0

Z = sum(rho[j, j] for j in range(nx + 1)) * dx
Z_an = sum(rho_an[j, j] for j in range(nx + 1)) * dx
pi_of_x = [rho[j, j] / Z for j in range(nx + 1)]
pi_of_x_an = [rho_an[j, j] / Z_an for j in range(nx + 1)]

pylab.plot(x, pi_of_x, c='red', linewidth=2)
pylab.plot(x, pi_of_x_an, c='blue', linewidth=0.5)
pylab.legend(["Harmonic", "Anharmonic"])
pylab.title('Cubic/Quartic Trotter Distribution $\pi(x)$ \
    \nfor Harmonic and Anharmonic Potentials for beta='+str(beta), fontsize = 18)
pylab.xlabel('$x$', fontsize = 30)
pylab.ylabel('$\pi(x)$', fontsize = 30)
pylab.show()
f = open('data_anharm_matrixsquaring_beta' + str(beta_tmp) + '.dat', 'w')
for j in range(nx + 1):
    f.write(str(x[j]) + ' ' + str(rho_an[j, j] / Z_an) + '\n')
f.close()
