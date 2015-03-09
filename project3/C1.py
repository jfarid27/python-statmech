from numpy import sum, cos, sin, tan, exp
from math import pi

# The angle of the circles off the axis
off_angle = tan(1/2)

angles_a = [0, (pi/2)-off_angle, (pi/2)+off_angle, pi, (3*pi/2)- off_angle, (3*pi/2)+off_angle]

#Just rotate the system by tan(1/2)
angles_b = map(lambda angle: angle + off_angle, angles_a)

#convert 15 degrees to radians
adjust = 15 * (pi/180)
#rotate by adjust angle
angles_c = map(lambda angle: angle + adjust, angles_a)

#group for processing
angles = [angles_a, angles_b, angles_c]

answers = [ sum(map(lambda p: exp(p*6*1j), angle_set))/6 for angle_set in angles]

print answers
