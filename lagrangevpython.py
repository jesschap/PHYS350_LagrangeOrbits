import pdb
import math
from vpython import *

# Simulation timestep, increasing this number
# will make the simulation run faster. Timestep
# currently set to a day. Note the sleep value is
# set to 0.1ms, which means that relative to a real second,
# the simulation will operate at timestep * 10000 seconds.
timestep = 24 * 3600

# Gravitational constant
G = 6.67428e-11

# Astronomical Unit
AU = 149.6e6 * 1000  # 149.6 million km, in meters.

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
E_MER = 0.205
E_VEN = 0.007
E_EAR = 0.017
E_MAR = 0.094
E_JUP = 0.049
E_SAT = 0.057
E_URA = 0.046
E_NEP = 0.011

# Angular Momentums of planets (kg*m^2/s)
L_MER =  8.95449 * 10**38
L_VEN = -1.84562 * 10**40
L_EAR =  2.66003 * 10**40
L_MAR =  3.51553 * 10**39
L_JUP =  1.92635 * 10**43
L_SAT =  7.82148 * 10**42
L_URA = -1.69313 * 10**42
L_NEP =  2.49140 * 10**42

# Distances of planets to sun (m)
# PROBLEM: THESE are al in AU so are mean distances. Change later
D_MER =  0.390 * AU
D_VEN =  0.723 * AU
D_EAR = -1.000 * AU
D_MAR = -1.524 * AU
D_JUP =  5.203 * AU
D_SAT =  9.539 * AU
D_URA =  19.18 * AU
D_NEP =  30.06 * AU

# Orbital tips of planets (converted from degrees to radians)
T_MER = 7.0 * (math.pi / 180)
T_VEN = 3.4 * (math.pi / 180)
T_EAR = 0.0 * (math.pi / 180)
T_MAR = 1.9 * (math.pi / 180)
T_JUP = 1.3 * (math.pi / 180)
T_SAT = 2.5 * (math.pi / 180)
T_URA = 0.8 * (math.pi / 180)
T_NEP = 1.8 * (math.pi / 180)


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


def compute_motion(planet):
    """
    planet: the planet to compute motion of

    This finds the new angle and distance of a planet from the sun
    """

    # The sun should never be passed as a parameter for this calculation
    if planet.name == 'Sun':
        raise ValueError()

    # Reduced mass constant
    u = (M_SUN * planet.mass) / (M_SUN + planet.mass)

    # Constant c value for the r(phi) function
    c = (planet.lz)**2 / (G * M_SUN * planet.mass * u)

    # Get the new planet angle
    d_angle = planet.lz / (u * planet.dist**2) * timestep
    planet.angle = planet.angle + d_angle

    # Get the new planet distance
    planet.dist = c / (1 + planet.ecc * math.cos(planet.angle))

    # Update the vpython simulation
    update_vmodel(planet)


def update_vmodel(planet):
    """
    planet: the planet to update vpython model of

    Update the vpython model of the planet by computing its x and y
    position from distance and angle.

    """
    x = planet.dist * math.cos(planet.angle)
    y = planet.dist * math.sin(planet.angle)

    # Update z component with tip formula
    z = math.sin(-planet.tip * math.cos(planet.angle)) * planet.dist

    # Assign the model position in x and y, assume it is at z = 0
    planet.model.pos = vector(x, y, z)


def loop(bodies):
    """
    bodies: a list of planets that will be modelled in the potential
            well of the sun

    Never returns; loops through the simulation, updating the
    positions of all the provided bodies.
    """
    while True:
        sleep(0.0001)
        for body in bodies:
            compute_motion(body)


def main():
#   pdb.set_trace()

    scene.title  = 'Planet Simulation'
    scene.width  = 1100
    scene.height = 700

    zplane = box(pos     = vector(0,0,0),
                 color   = color.white,
                 length  = D_NEP,
                 height  = D_NEP,
                 width   = R_SUN / 10,
                 opacity = 0.5)

    # Initialize all of the vpython models
    sunmodel     = sphere(pos = vector(0,0,0),
                          radius = R_SUN * SUNSCALE,
                          color = color.yellow)
    mercurymodel = sphere(pos = vector(D_MER,0,0),
                          radius = R_MER * SMALLBODYSCALE,
                          color = vec(1,1,1))
    venusmodel   = sphere(pos = vector(D_VEN,0,0),
                          radius = R_VEN * SMALLBODYSCALE,
                          color = color.orange)
    earthmodel   = sphere(pos = vector(D_EAR,0,0),
                          radius = R_EAR * SMALLBODYSCALE,
                          color = color.blue)
    marsmodel    = sphere(pos = vector(D_MAR,0,0),
                          radius = R_MAR * SMALLBODYSCALE,
                          color = color.red)
    jupitermodel = sphere(pos = vector(D_JUP,0,0),
                          radius = R_JUP * LARGEBODYSCALE,
                          color = color.green)
    saturnmodel  = sphere(pos = vector(D_SAT,0,0),
                          radius = R_SAT * LARGEBODYSCALE,
                          color = vec(1,0,1))
    uranusmodel  = sphere(pos = vector(D_URA,0,0),
                          radius = R_URA * LARGEBODYSCALE,
                          color = color.white)
    neptunemodel = sphere(pos = vector(D_NEP,0,0),
                          radius = R_NEP * LARGEBODYSCALE,
                          color = color.blue)

    # Initialize all of the planet objects
    sun     = Body('Sun', sunmodel, M_SUN, 0, 0, 0, 0, 0)
    mercury = Body('Mercury', mercurymodel, M_MER, 0, D_MER, L_MER, E_MER, T_MER)
    venus   = Body('Venus', venusmodel, M_VEN, 0, D_VEN, L_VEN, E_VEN, T_VEN)
    earth   = Body('Earth', earthmodel, M_EAR, 0, D_EAR, L_EAR, E_EAR, T_EAR)
    mars    = Body('Mars', marsmodel, M_MAR, 0, D_MAR, L_MAR, E_MAR, T_MAR)
    jupiter = Body('Jupiter', jupitermodel, M_JUP, 0, D_JUP, L_JUP, E_JUP, T_JUP)
    saturn  = Body('Saturn', saturnmodel, M_SAT, 0, D_SAT, L_SAT, E_SAT, T_SAT)
    uranus  = Body('Uranus', uranusmodel, M_URA, 0, D_URA, L_URA, E_URA, T_URA)
    neptune = Body('Neptune', neptunemodel, M_NEP, 0, D_NEP, L_NEP, E_NEP, T_NEP)

    loop([mercury, venus, earth, mars, jupiter, saturn, uranus, neptune])


if __name__ == '__main__':
    main()
