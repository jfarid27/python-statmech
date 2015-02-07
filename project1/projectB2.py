import random

deltas = [ 0.062, 0.125, 0.25, 0.5, 1, 2, 4 ] 

def computeAcceptance(delta, n_trials=int(2e6)):

    x, y = 1.0, 1.0
    accepted = 0

    for i in range(n_trials):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
            accepted += 1
    
    return delta, accepted / float(n_trials)

if __name__ == '__main__':

    acceptances = []

    for delta in deltas:
        print computeAcceptance(delta)
