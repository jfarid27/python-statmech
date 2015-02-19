import random, math

n_trials = 400000
obs = 0
sq_obs = 0
for iter in range(n_trials):
    x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
    if x**2 + y**2 < 1.0:
        obs += 4.0
        sq_obs += 16

avg_obs = obs /float(n_trials)
avg_obs_sq = sq_obs / float(n_trials)
var =  avg_obs_sq - avg_obs**2 
sd = math.sqrt(var)
print sd
