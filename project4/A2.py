import random

x, y = 0.0, 0.0
delta = 0.1
n_trials = 1000000
n_hits = 0.0
for i in range(n_trials):
    del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
    #forcing sampling in the unit circle here
    if abs(x+del_x)**2 + abs(y+del_y)**2 < 1.0: 
        x, y = x + del_x, y + del_y
    z = random.uniform(-1.0,1)
    if x**2 + y**2 + z**2 < 1.0: n_hits += 1

#not using trials here due to the forcing of sampling in the unit circle
print 2.0 * float(n_hits) / float(n_trials)
