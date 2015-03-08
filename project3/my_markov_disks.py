import random, math, os, pylab
from itertools import product

def dist(x,y):
    d_x = abs(x[0] - y[0]) % 1.0
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(x[1] - y[1]) % 1.0
    d_y = min(d_y, 1.0 - d_y)
    return  math.sqrt(d_x**2 + d_y**2)

def extend_periodic(disk_configurations, periods=[1,1]):
    """
    Function to extend disks into periodic boundary conditions.
    Default periods array is an array parameter with values
    specifying the period in each dimension of the coordinate system.
    Hence the default behavior extends into the first dimension and
    second dimensions with period 1. It is assumed each disk in disk
    in disk configuration is an array with length matching periods
    array length.
    """
    configuration = []
    
    #making all the positions in each dim
    enumerations = [[-range, 0, range] for range in periods]
    #taking the cartesian product of all possible positions
    zero_points = [list(p) for p in product(*enumerations)]
    for disk in disk_configurations:
        #create a new position set representing the disk
        #in periodic system
        disk_positions = list(zero_points)
        #add disk position to periodic system set
        for ind, point in enumerate(disk_positions):
            for dim, pos in enumerate(point):
                disk_positions[ind][dim] += disk[dim]
        #add all new disk positions to configuration
        for new_disk in disk_positions:
            configuration.append(new_disk)
    return configuration

def show_conf(L, sigma, title, fname):
    """
    Function to plot given configuration of circles.
    """
    pylab.axes()
    for [x, y] in L:
        cir = pylab.Circle((x, y), radius=sigma,  fc='r')
        pylab.gca().add_patch(cir)
    pylab.axis('scaled')
    pylab.title(title)
    pylab.axis([0.0, 1.0, 0.0, 1.0])
    pylab.savefig(fname)
    pylab.show()
    pylab.close()

def initial_conf(default=False, filename="disk_configuration.txt", N=4,sigma=.1):
    
    n_sqrt = int(math.sqrt(N))
    delxy = sigma
    two_delxy = sigma * 2
    if not (default):
        return [[delxy + i * two_delxy, delxy + j * two_delxy] for i in range(n_sqrt) for j in range(n_sqrt)]
    else:
        L = []
        if os.path.isfile(filename):
            f = open(filename, 'r')
            for line in f:
                a, b = line.split()
                L.append([float(a), float(b)])
            f.close()
            print 'starting from file', filename
            return L
        else:
            L = [[delxy + i * two_delxy, delxy + j * two_delxy] for i in range(n_sqrt) for j in range(n_sqrt)]
            print 'starting from a new random configuration'
        f = open(filename, 'w')
        for a in L:
            f.write(str(a[0]) + ' ' + str(a[1]) + '\n') 
        return L

def run_chain(eta=.1, N=4, plot=False, png="sim_output", inconf=False):
    """
    Function to run the markov chain simulation with optional density
    eta for four disks. If plot parameter is true, function will
    output plot to specified file name.
    """
    #Calculate radius from specified density
    sigma = math.sqrt( eta / (math.pi * N) )
    sigma_sq = sigma ** 2
    delta = 0.1
    n_steps = 10
    filename = 'disk_configuration_N%i_eta%.2f.txt' % (N, eta)
    L = initial_conf(inconf, filename, N, sigma)
    for steps in range(n_steps):
        #Choose a disk
        a = random.choice(L)
        #update the disks position
        b = [a[0] + random.uniform(-delta, delta), a[1] + random.uniform(-delta, delta)]
        #Expand the disks using periodic boundary conditions
        expansions = extend_periodic(L)
        #Do minimum distance check
        min_dist = min([dist(b,c) for c in expansions if c != a])
        if not (min_dist < sigma * 2):
            #fold disks back into interval
            new_a = [ pos % 1 if pos > 1 or pos < 0 else pos for pos in b]
            a[:] = new_a
    if (plot == True):
        exps = extend_periodic(L,[1,1])
        show_conf(exps, sigma, "Simulation", png + ".png")
    return L

if __name__ == "__main__":
    conf = run_chain(eta=.7853, N=64, plot=True, inconf=True)
