import random, math

def markov_compression(d, delta=0.1, n_trials=100000):
    if d == 2:
        return math.pi/2
    #The calculation for the Q(d) value
    conf = [0.0 for x in range(0,d)]
    n_hits = 0.0
    n_samples = 0.0
    for i in range(n_trials):
        k = random.randint(0, d - 1)
        x_old_k = conf[k]
        x_new_k = x_old_k + random.uniform(-delta, delta)
        if abs(x_new_k) < 1:
            #forcing sampling in the unit circle
            new_rad = sum([x**2 for x in conf]) - x_old_k**2 + x_new_k**2 
            if new_rad < 1: 
                n_samples += 1
                alpha = random.uniform(-1.0,1.0)
                if new_rad + alpha**2 < 1:
                    n_hits +=1
            #update configuration with new sample
            conf[k] = x_new_k
    return 2 * n_hits / float(n_samples)

def product_vol(d):
    return 2 * reduce(lambda a, b: a*b, [markov_compression(dim) for dim in range(2, d+1)])


if __name__ == "__main__":
    print product_vol(3)
