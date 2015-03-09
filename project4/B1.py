import random, math

def markov_hypersphere(d, delta=0.1, n_trials=1000000):
    if d ==1:
        return 2
    conf = [0.0 for x in range(0,d)]
    n_hits = 0.0
    for i in range(n_trials):
        k = random.randint(0, d - 1)
        x_old_k = conf[k]
        x_new_k = x_old_k + random.uniform(-delta, delta)
        if abs(x_new_k) < 1:
            #forcing sampling in the unit circle
            if sum([x**2 for x in conf]) - x_old_k**2 + x_new_k**2  < 1.0: 
                n_hits += 1
            #update configuration with new sample
            conf[k] = x_new_k
    return 4 * n_hits / float(n_trials)
