import random, math

class StepControlledMCMC():

    def prob(self, x):
        s1 = math.exp(-(x + 1.2) ** 2 / 0.72)
        s2 = math.exp(-(x - 1.5) ** 2 / 0.08)
        return (s1 + 2.0 * s2) / math.sqrt(2.0 * math.pi)

    def controlDelta(self, prev_delta, acc_rate):
        """Returns new delta to control the temporary acceptance rate.
        """
        if (acc_rate > .6):
            return prev_delta * 1.1
        elif (acc_rate < .4):
            return prev_delta / 1.1
        else:
            return prev_delta

    def simulate(self, delta=10.0, nsteps=10000, hist=False):
        """Returns a tuple of a simulation's acceptance rate and averages
           with specified delta and nsteps. If hist named argument is true,
           the function also returns a list of all the x values.
        """
        acc_tot = 0
        acc_temp = 0
        temp_size = 100
        x = 0.0
        x_av = 0.0
        x_vals = []
        for step in xrange(nsteps):
            if (step % temp_size == 0):
                #Modification to control the delta
                delta = self.controlDelta(delta, acc_temp/float(temp_size))
                acc_temp = 0
            xnew = x + random.uniform(-delta, delta)
            if random.uniform(0.0, 1.0) < self.prob(xnew) / self.prob(x):
                x = xnew
                acc_tot += 1
                acc_temp += 1
                if hist: x_vals.append(x)
            x_av += x
        if (hist):
            return acc_tot / float(nsteps), x_av / float(nsteps), x_vals
        return acc_tot / float(nsteps), x_av / float(nsteps)
        

if __name__ == "__main__":
    max_order_mag = 6
    deltas=[10 ** x for x in range(-max_order_mag, max_order_mag, 1)]
    nsteps = 1000000
    inst = StepControlledMCMC()
    for delta in deltas:
        acc_rate, x_av = inst.simulate(delta, nsteps)
        print "-----------------------"
        print "delta = %f" % delta
        print "acc_rate = %f" % acc_rate
        print "x_av = %f" % x_av 
