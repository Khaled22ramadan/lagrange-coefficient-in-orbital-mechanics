import numpy as np
import math 


# at first i entred the initial value of the position and velocity vector 

# Initial position vector
x = [8182.4, -6865.9, 0] 
R_initial = np.array(x)         
x_mag = np.linalg.norm(R_initial)
print ("the position vector X = : %s" % R_initial)           #print (R_initial)
print ("Magnitude of position vector = %f" % x_mag)          #print (V_initial)
print ("-----------------------------------")

# Initial velocity vector
v = [0.47572, 8.8116, 0]
V_initial = np.array(v)         
v_mag = np.linalg.norm(V_initial)
print ("the velocity vector V = : %s" % V_initial)           #print (V_initial)
print ("Magnitue of velocity vector = : %f" % v_mag)
print ("-----------------------------------")

#true anomly
angle  = 120                    
angle_in_radians = math.radians(120)                        # Convert angle from degree to radians
print ("the true anomly = %d " % angle)
# sin & cos values 
sin_value = np.sin(angle_in_radians)
cos_value = np.cos(angle_in_radians)
print ("cos value = %f " %cos_value)
print ("sin value = %f " %sin_value)
print ("-----------------------------------")

# the result of the dot product  for X&V initial and it's magnitude
rovodot = dot_product = np.dot(R_initial, V_initial)  #dot process 
print ("rovo is the result for dot produdt for r&v initial = %f" %rovodot)
print ("-----------------------------------")

# the result of the cross product  for X&V initial and it's magnitude (h)
h = cross_product = np.cross(R_initial, V_initial)    # cross porcess 
h_mag = np.linalg.norm(cross_product)  
print ("the result of cross product for r&v = %s" % h)
print ("Magnitude of h = %f" %h_mag)
print ("-----------------------------------")

#initial radial velocity vr0
Vro = (rovodot/ x_mag)
print ("Initial radial velocity = %f" %Vro)
print ("-----------------------------------")

#The final distance r ==> (h^2 /M) * (1/ 1+(h^2/Mro-1)*cos(ang) - (h*Vro/M)*sin(ang))

hspec = 75366 # specific angular momentum const
M = 398600    # const 
r = (hspec**2 / M) * (1 / (1 + (((hspec**2 / (M * x_mag)) - 1) * cos_value) - ((hspec * Vro * sin_value) / M )))
#r = (hspec**2/ M) * (1)/ ( 1+ ( hspec**2 / (M * x_mag) _1) *cos_value) _ ((hspec *Vro)/M)*sin_value )
print ("the final distance r = %f" %r)
print ("-----------------------------------")

# now i will calc lagrange cooeficient f ,g , f. & g.
print ("now i will calculate lagrange cooeficients f ,g , fdot & gdot")
#f
f = 1 - (((M * r) / hspec**2 ) * (1 - cos_value))
print ("f = %f" %f)

#g
g = (r * x_mag * sin_value) / hspec
print ("g = %f" %g)

# fdot
ff = ((M * (1 - cos_value)) / (hspec * sin_value)) * (((M * (1 - cos_value))/ hspec**2 )- (1/x_mag) - (1/ r))  
print ("fdot = %f" %ff)

# gdot
gg = 1 - ((M * x_mag) * (1 - cos_value) / hspec**2 )
print ("gdot = %f "%gg)
print ("-----------------------------------")

# now we can calculate the value of the position and velocity for space craft at any time in orbit  using initial values for them

# find the final position vector (R_new)
f = np.array(f)    # (f) is converted to arrary to multibly it in initial position vector 
g = np.array(g)    # (g) is converted to arrary to multibly it in initial velocity vector 
R_new = (f * x) + (g * v)
print ("the final position vector = %s" % R_new)

# find the final velocity vector (V_new)
ff = np.array(ff)   # ff is the (fdot) and i converted to arrary to multibly it in initial position vector 
gg = np.array(gg)   # ff is the (gdot) and i converted to arrary to multibly it in initial velocity vector 
V_new = (ff * x) + (gg * v)
print ("the final velocity vector = %s" % V_new)

print ("*********************")




