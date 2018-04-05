# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 16:43:29 2018

@author: Aaron Janz
       : Jessica Chapman
       : Declan Johnston
       : Theophilus Ko
"""
import pdb
import math
import random as rnd
from vpython import *

# Simulation timestep, increasing this number
# will make the simulation run faster. Timestep
# currently set to a day. Note the sleep value is
# set to 0.1ms, which means that relative to a real second,
# the simulation will operate at timestep * 10000 seconds.
timestep = 24*3600*2

collisions = 0
# Gravitational constant
G = 6.67428e-11

# Astronomical unit (converted to meters)
AU = 149.6e6 * 1000

# Scale factor for planets size so that they're visible
# relative to the large distances between them.
SUNSCALE = 50          # For use by Sun, it's too large
LARGEBODYSCALE = 500   # For use by the giants, they're too large
SMALLBODYSCALE = 2000  # For use by other planets

# Masses of all the planets (kg)
M_SUN = 1.989 * 10**30
M_MER = 3.3   * 10**23
M_VEN = 4.87  * 10**24
M_EAR = 5.97  * 10**24
M_MAR = 6.42  * 10**23
M_JUP = 1.898 * 10**27
M_SAT = 5.68  * 10**26
M_URA = 8.68  * 10**25
M_NEP = 1.02  * 10**26

# Radiuses of all the bodies (m)
R_SUN = 696000 * 10**3
R_MER = 2448.5 * 10**3
R_VEN = 6052   * 10**3
R_EAR = 6378   * 10**3
R_MAR = 3396   * 10**3
R_JUP = 71492  * 10**3
R_SAT = 60268  * 10**3
R_URA = 25559  * 10**3
R_NEP = 24764  * 10**3

# Eccentricities of planets (unitless)
E_SUN = 0.000
E_MER = 0.205
E_VEN = 0.007
E_EAR = 0.017
E_MAR = 0.094
E_JUP = 0.049
E_SAT = 0.057
E_URA = 0.046
E_NEP = 0.011

# Angular Momentums of planets (kg*m^2/s)
L_SUN =  5.95498 * 10**41
L_MER =  8.95449 * 10**38
L_VEN = -1.84562 * 10**40
L_EAR =  2.66003 * 10**40
L_MAR =  3.51553 * 10**39
L_JUP =  1.92635 * 10**43
L_SAT =  7.82148 * 10**42
L_URA = -1.69313 * 10**42
L_NEP =  2.49140 * 10**42

# Distances of planets to sun (m)
D_SUN =  0.7420 * 10**9
D_MER =  57.900 * 10**9
D_VEN =  108.20 * 10**9
D_EAR = -149.60 * 10**9
D_MAR = -227.90 * 10**9
D_JUP =  778.60 * 10**9
D_SAT =  1433.5 * 10**9
D_URA =  2872.5 * 10**9
D_NEP =  4495.1 * 10**9

# Initial Phase Shifts of Planets (conveted from degrees to radians)
P_MER = 140.0 * (math.pi / 180)
P_VEN =  45.0 * (math.pi / 180)
P_EAR = 180.0 * (math.pi / 180)
P_MAR = 230.0 * (math.pi / 180)
P_JUP = 225.0 * (math.pi / 180)
P_SAT = 280.0 * (math.pi / 180)
P_URA = 340.0 * (math.pi / 180)
P_NEP =  30.0 * (math.pi / 180)

# Orbital tips of planets (converted from degrees to radians)
T_MER = 7.0 * (math.pi / 180)
T_VEN = 3.4 * (math.pi / 180)
T_EAR = 0.0 * (math.pi / 180)
T_MAR = 1.9 * (math.pi / 180)
T_JUP = 1.3 * (math.pi / 180)
T_SAT = 2.5 * (math.pi / 180)
T_URA = 0.8 * (math.pi / 180)
T_NEP = 1.8 * (math.pi / 180)

# Asteroid Data
M_AST  = 10**16
M_AST1 = 2.2 * 10**14
M_AST2 = 10  * 10**14

R_AST1 = 10000 * 10**3
R_AST2 = 20000 * 10**3

E_AST1 = 1.0
E_AST2 = 0.6

D_AST  = 40.0 * AU
D_AST1 = 12.0 * AU
D_AST2 = 10.0 * AU

L_AST1 = M_AST1*D_AST1*5560
L_AST2 = M_AST2*D_AST1*8000

class Body:
    """
    name  : name of planet
    model : vpython model of planet
    mass  : mass (kg)
    angle : planet angle from reference (radians)
    dist  : distance to sun (m)
    lz    : z angular momentum (kgm^2/s)
    ecc   : eccentricity (unitless)
    """

    def __init__(self, name, model, mass, angle, dist, lz, ecc, tip):
        self.name  = name
        self.model = model
        self.mass  = mass
        self.angle = angle
        self.dist  = dist
        self.lz    = lz
        self.ecc   = ecc
        self.tip   = tip


class Asteroid:
    """
    name    : name of planet
    model   : vpython model of planet
    mass    : mass (kg)
    angle   : planet angle from reference (radians)
    vx, vy  : x, y velocities in m/s
    px, py  : x,y positions in m
    ivx, ivy: initial x, y velocities in m/s
    ipx, ipy: initial x,y positions in m
    collided: true if the asteroid has collided with earth
    """

    def __init__(self, name, mass, angle, vx, vy, px, py):
        self.name  = name
        self.mass  = mass
        self.angle = angle
        self.vx    = vx
        self.vy    = vy
        self.px    = px
        self.py    = py
        self.ivx   = vx
        self.ivy   = vy
        self.ipx   = px
        self.ipy   = py
        self.model = asteroidmodel1 = sphere(pos    = vector(px,py,0),
                                radius = R_AST1 * LARGEBODYSCALE,
                                color  = color.cyan, make_trail=True, retain = 100)
        self.collided = False

    def attraction(self, other,asteroids):
        """(Body): (fx, fy)

        Returns the force exerted upon this body by the other body.
        """
        # Report an error if the other object is the same as this one.
        if self is other:
            raise ValueError("Attraction of object %r to itself requested"
                             % self.name)

        # Compute the distance of the other body.
        sx, sy = self.px, self.py
        #print("Asteroid x pos: %f" %(sx/AU))
        #print("Asteroid y pos: %f" %(sy/AU))
        ox, oy = other.dist*math.cos(other.angle), other.dist*math.sin(other.angle)
        dx = (ox-sx)
        dy = (oy-sy)
        d = math.sqrt(dx**2 + dy**2)


        # Report an error if the distance is zero cuz u
        # get a ZeroDivisionError exception further down.
        if d <= 2*R_EAR and other.name == 'Earth':
            global collisions
            collisions += 1
            self.model.color = color.red
            asteroids.remove(self)

        # If its within a certain radius of the earth, it's collided. Keep
        # track of collisions so only record one collision per asteroid.
        # Print out its initial conditions to the screen
        if d <= AU*0.5 and other.name == 'Earth' and self.collided == False:
            global collisions
            collisions += 1
            self.collided = True
            print("\nInitial position: " + str(self.ipx/AU) + ", " + str(self.ipy/AU))
            print("Initial velocity: " + str(self.ivx) + ", " + str(self.ivy))
            f = open("collision_log_2R_EARTH_Jupiter_True.txt","a")
            f.write("\nInitial position in AU (x,y): " + str(self.ipx/AU) + ", " + str(self.ipy/AU))
            f.write("\nInitial velocity in AU (x,y): " + str(self.ivx) + ", " + str(self.ivy))
            f.write("\n")
            f.close()



        # Compute the force of attraction
        f = G * self.mass * other.mass / (d**2)

        # Compute the direction of the force.
        theta = math.atan2(dy, dx)
        #print("theta: %f" %theta)
        fx = math.cos(theta) * f
        fy = math.sin(theta) * f

        return fx, fy



def compute_forces(ast,bodies,asteroids):

    total_fx = total_fy = 0.0
    for body in bodies:
        # Add up all of the forces exerted on 'body'.
        fx, fy = ast.attraction(body,asteroids)
        total_fx += fx
        total_fy += fy

    # Update velocities based upon on the force.
    ast.vx += total_fx / ast.mass * timestep
    ast.vy += total_fy / ast.mass * timestep

    # Update positions
    ast.px += ast.vx * timestep
    ast.py += ast.vy * timestep

    #Updates positions of asteroids in simulation
    ast.model.pos = vector(ast.px,ast.py,0)

def compute_motion(body):
    """
    body: the body to compute the motion of

    This finds the new angle and distance of a body from the center of
    the system's mass
    """

    # Compute this assuming the sun is the main mass
    M_OBJ = M_SUN

    # If computing motion of sun, use Jupiter as its mass
    if body.name == 'Sun':
        M_OBJ = M_JUP

    # Reduced mass constant
    u = (M_OBJ * body.mass) / (M_OBJ + body.mass)

    # Constant c value for the r(phi) function
    c = (body.lz)**2 / (G * M_OBJ * body.mass * u)

    # Get the new body angle
    d_angle = body.lz / (u * body.dist**2) * timestep
    body.angle = body.angle + d_angle

    # Get the new body distance
    body.dist = c / (1 + body.ecc * math.cos(body.angle))

    # Update the vpython simulation
    update_vmodel(body)


def update_vmodel(body):
    """
    body: the body to update vpython model of

    Update the vpython model of the body by computing its x and y
    position from distance and angle.

    """
    x = body.dist * math.cos(body.angle)
    y = body.dist * math.sin(body.angle)

    # Update z component with tip formula
    z = math.sin(-body.tip * math.cos(body.angle)) * body.dist

    # Assign the model position in x and y, assume it is at z = 0
    body.model.pos = vector(x, y, z)

def loop(bodies, asteroids):
    """
    bodies: a list of bodies that will be modelled in the potential
            well of the sun

    Never returns; loops through the simulation, updating the
    positions of all the provided bodies.
    """
    count = 0
    
    while True:
        if count % 1 == 0:
            # Print the years lapsed with commas to separate 3's of digits.
            scene.title  = 'Planet Simulation: ' + str("{:,d}".format(int(timestep*count/(3.152*10**7))) + 
                                                       ' years' + '   Collisions: ' +str(collisions) +
                                                       '   Number of Local Asteroids: ' +str(len(asteroids)))
#            print(str("{:,d}".format(int(timestep*count/365))) + ' years')
        if len(asteroids) < 30:
            newAsteroid = Asteroid('Asteroid'+str(count),
                                       M_AST*float(rnd.randint(0,10)+rnd.uniform(-1,1)), 0, 10000*float(rnd.uniform(-1,1)),10000*float(rnd.uniform(-1,1)),
                                       (30*AU+rnd.randint(-3,3)*AU)*(-1)**rnd.randrange(2),(30*AU+rnd.randint(-3,3)*AU)*(-1)**rnd.randrange(2))
            asteroids.append(newAsteroid)
        asteroidDummy = list(asteroids)
        #sleep(0.0001)
        for body in bodies:
            compute_motion(body)
        for ast in asteroidDummy:
            if abs(ast.px) > 100*AU or abs(ast.py) > 100*AU:
                ast.model.retain = 1
                asteroids.remove(ast)

            compute_forces(ast,bodies,asteroids)
        count += 1
def main():
#   pdb.set_trace()

    scene.title  = 'Planet Simulation'
    scene.width  = 1100
    scene.height = 700
    #text(text='Years', align = 'center')
    # Inserting this ugly plane will help us see orientation better. We
    # can likely read input from user to toggle this plane on/off by setting
    # it's opacity level
    zplane = box(pos     = vector(0,0,0),
                 color   = color.white,
                 length  = 2 * D_NEP,
                 height  = 2 * D_NEP,
                 width   = R_SUN / 10,
                 opacity = 0.3)

    # Initialize all of the vpython models
    sunmodel     = sphere(pos    = vector(D_SUN,0,0),
                          radius = R_SUN * SUNSCALE,
                          texture={'file':textures.flower},
                          color  = color.yellow)
    mercurymodel = sphere(pos    = vector(D_MER,0,0),
                          radius = R_MER * SMALLBODYSCALE,
                          color  = vec(1,1,1))
    venusmodel   = sphere(pos    = vector(D_VEN,0,0),
                          radius = R_VEN * SMALLBODYSCALE,
                          texture={'file':textures.stucco},
                          color  = color.orange)
    earthmodel   = sphere(pos    = vector(D_EAR,0,0),
                          radius = R_EAR * SMALLBODYSCALE,
                          texture={'file':textures.earth},
                          make_trail = True,
                          retain = 5)
    marsmodel    = sphere(pos    = vector(D_MAR,0,0),
                          radius = R_MAR * SMALLBODYSCALE,
                          texture={'file':textures.stucco},
                          color  = color.red)
    jupitermodel = sphere(pos    = vector(D_JUP,0,0),
                          radius = R_JUP * LARGEBODYSCALE,
                          texture={'file':textures.stucco},
                          color  = color.green)
    saturnmodel  = sphere(pos    = vector(D_SAT,0,0),
                          radius = R_SAT * LARGEBODYSCALE,
                          texture={'file':textures.stucco},
                          color  = vec(1,0,1))
    uranusmodel  = sphere(pos    = vector(D_URA,0,0),
                          radius = R_URA * LARGEBODYSCALE,
                          texture={'file':textures.stucco},
                          color  = color.white)
    neptunemodel = sphere(pos    = vector(D_NEP,0,0),
                          radius = R_NEP * LARGEBODYSCALE,
                          texture={'file':textures.stucco},
                          color  = color.blue)

    # Initialize all of the planet objects
    sun     = Body('Sun',     sunmodel,     M_SUN, 0,     D_SUN, L_SUN, E_SUN, 0)
    mercury = Body('Mercury', mercurymodel, M_MER, P_MER, D_MER, L_MER, E_MER, T_MER)
    venus   = Body('Venus',   venusmodel,   M_VEN, P_VEN, D_VEN, L_VEN, E_VEN, T_VEN)
    earth   = Body('Earth',   earthmodel,   M_EAR, P_EAR, D_EAR, L_EAR, E_EAR, T_EAR)
    mars    = Body('Mars',    marsmodel,    M_MAR, P_MAR, D_MAR, L_MAR, E_MAR, T_MAR)
    jupiter = Body('Jupiter', jupitermodel, M_JUP, P_JUP, D_JUP, L_JUP, E_JUP, T_JUP)
    saturn  = Body('Saturn',  saturnmodel,  M_SAT, P_SAT, D_SAT, L_SAT, E_SAT, T_SAT)
    uranus  = Body('Uranus',  uranusmodel,  M_URA, P_URA, D_URA, L_URA, E_URA, T_URA)
    neptune = Body('Neptune', neptunemodel, M_NEP, P_NEP, D_NEP, L_NEP, E_NEP, T_NEP)
    #(name, model, mass, angle,vx, vy, px, py)
    #10000.0*float(rnd.uniform(-1,1))
    asteroid1 = Asteroid('Asteroid1', M_AST1, 0, -5000.0,0, 10*AU, 10*AU)
    asteroid2 = Asteroid('Asteroid2', M_AST2, 0, 3000, 0, -10,-10*AU)
    #First quadrant-> negative vx and vy
    asteroid3 = Asteroid('Asteroid3', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 7000.0*rnd.uniform(-1,-0.5),7000.0*rnd.uniform(-1,-0.5), 20*AU,30*AU)
    asteroid4 = Asteroid('Asteroid4', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 17000.0*rnd.uniform(-1,-0.5),17000.0*rnd.uniform(-1,-0.5), 30*AU,20*AU)
    #Second Quadrant -> positive vx and negative vy
    asteroid5 = Asteroid('Asteroid5', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 7000.0*rnd.uniform(0.5,1),7000.0*rnd.uniform(-1,-0.5), -20*AU,30*AU)
    asteroid6 = Asteroid('Asteroid6', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 17000.0*rnd.uniform(0.5,1),17000.0*rnd.uniform(-1,-0.5), -30*AU,20*AU)
    #Third Quadrant -> positive vx and vy
    asteroid7 = Asteroid('Asteroid7', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 7000.0*rnd.uniform(0.5,1),7000.0*rnd.uniform(0.5,1), -25*AU,-30*AU)
    asteroid8 = Asteroid('Asteroid8', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 17000.0*rnd.uniform(0.5,1),17000.0*rnd.uniform(0.5,1), -30*AU,-25*AU)
    #4th Quadrant -> negative vx and positive vy
    asteroid9 = Asteroid('Asteroid9', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)), 0, 7000.0*rnd.uniform(-1,-0.5),7000.0*rnd.uniform(0.5,1), 27*AU,-30*AU)
    asteroid10= Asteroid('Asteroid10', M_AST*(rnd.randint(1,100)+rnd.uniform(-1,1)),0, 17000.0*rnd.uniform(-1,-0.5),17000.0*rnd.uniform(0.5,1), 30*AU,-27*AU) 

    loop([mercury, venus,sun, earth, mars,jupiter,saturn, uranus, neptune], [])

if __name__ == '__main__':
    main()




