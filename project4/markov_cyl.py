import random

def markov_cyl(n_trials=100000, delta=.75, debug=False):
  x, y, z = 0.0, 0.0, 0.0
  n_hits = 0
  accepted = 0
  for i in range(n_trials):
      del_x = random.uniform(-delta, delta) 
      del_y = random.uniform(-delta, delta)
      del_z = random.uniform(-delta, delta)
      if abs(x+del_x) < 1 and abs(y+del_y) < 1.0 and abs(z+del_z) < 1.0: 
          accepted += 1
          x, y, z = x + del_x, y + del_y, z + del_z
      if x**2 + y**2 < 1.0:
          n_hits += 1
  ans = 8 * float(n_hits) / n_trials
  if debug:
    print accepted / float(n_trials)
    print ans
  return ans

if __name__ == "__main__":
  markov_cyl(debug=True)
