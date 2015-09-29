import math, random, time

class SimAnnealing():

    def V(self, x, y):
        """Returns computed potential function.

        """
        pot  = -4.0 * x ** 2 - x ** 3 + 4.0 * x ** 4
        pot += -4.0 * y ** 2 - y ** 3 + 4.0 * y ** 4
        return pot

    def nearMinimum(self, x, y, tolerance = 0.1):
        """Returns true if x and y are near expected minimum within
           specified tolerance
        """
        xmin, ymin = 0.807044513157, 0.807044513157
        if ( (x - xmin) ** 2 + (y - ymin) ** 2 < tolerance):
            return True
        return False

    def simulate(self, gamma=0.4, n_runs=100) :
        """Returns simulation success rate and final temperature 
           results for specified number of runs and gamma
        """
        n_success = 0
        wall_clock = 0
        for run in range(n_runs):
            T = 4.0
            x, y = 0.0, 0.0
            delta = 0.1
            step = 0
            acc = 0
            start = time.time()
            while T > 0.00001:
                step += 1
                if step == 100:
                    T *= (1.0 - gamma)
                    if acc < 30:
                       delta /= 1.2
                    elif acc > 70:
                       delta *= 1.2
                    step = 0
                    acc = 0
                xnew = x + random.uniform(-delta, delta)
                ynew = y + random.uniform(-delta, delta)
                if abs(xnew) < 1.0 and abs(ynew) < 1.0 and \
                   random.uniform(0.0, 1.0) < math.exp(- (self.V(xnew, ynew) - self.V(x, y)) / T):
                    x = xnew
                    y = ynew
                    acc += 1
            if self.nearMinimum(x, y):
                n_success += 1
            wall_clock += time.time() - start 
        return  T, n_success / float(n_runs), wall_clock / float(n_runs)


if __name__ == "__main__":
    gammas=[2 ** x for x in range(-7, -1, 1)]
    nruns = 100
    inst = SimAnnealing()
    for gamma in gammas:
        final_temp, success_rate, avg_time = \
            inst.simulate(gamma, nruns)
        print "-----------------------"
        print "gamma = %f" % gamma
        print "success_rate = %f" % success_rate
        print "final_temp = %f" % final_temp
        print "Avg wall clock time = %f" % avg_time
