import random

x, y, z = 0.0, 0.0, 0.0
delta = 0.1
n_trials = 1000000
n_hits = 0.0
n_samples = 0.0
for i in range(n_trials):
    del_x = random.uniform(-delta, delta)
    del_y = random.uniform(-delta, delta)
    del_z = random.uniform(-delta, delta)
    if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0 and abs(z + del_z) < 1.0:
        #forcing sampling in the unit circle here
        if x**2 + y**2 + z**2 < 1.0: 
            n_samples += 1
            alpha = random.uniform(-1.0,1.0)
            if x**2 + y**2 + z**2 + alpha**2 < 1.0: n_hits += 1
        x, y, z = x + del_x, y + del_y, z + del_z

#not using trials here due to the forcing of sampling in the unit circle
print n_hits / float(n_samples)
