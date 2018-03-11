# pylama:ignore=W391,E231,E251,E303,E221,E272

import math
from vpython import *

# Gravitational constant
G = 6.67428e-11

# Astronomical Unit
AU = 149.6e6 * 1000  # 149.6 million km, in meters.

# Scale factor for planets size so that they're visible relative to the large
# distances between them.
SUNSCALE = 70  # For use by Sun, it's too large
LARGEBODYSCALE = 500  # For use by the giants, they're too large
SMALLBODYSCALE = 2000 # For use by other planets

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

class Body():
    """
    mass : mass in kg
    vx, vy: x, y velocities in m/s
    px, py: x, y positions in m
    """

    name = None
    mass = None
    model = None
    vx = vy = 0.0
    px = py = 0.0

    # This function computes the force between one body and the other, and then
    # returns them in the order fx, fy
    def attraction_force(self, other):
        """(Body): (fx, fy)

        Returns the force exerted upon this body by the other body.
        """
        # Report an error if the other object is the same as this one.
        if self is other:
            raise ValueError("Attraction of object %r to itself requested"
                             % self.name)

        # Compute the distance of the other body.
        sx, sy = self.px, self.py
        ox, oy = other.px, other.py
        dx = (ox-sx)
        dy = (oy-sy)
        d = math.sqrt(dx**2 + dy**2)

        # Report an error if the distance is zero cuz u
        # get a ZeroDivisionError exception further down.
        if d == 0:
            raise ValueError("Collision between objects %r and %r"
                             % (self.name, other.name))

        # Compute the force of attraction
        f = G * self.mass * other.mass / (d**2)

        # Compute the direction of the force.
        theta = math.atan2(dy, dx)
        fx = math.cos(theta) * f
        fy = math.sin(theta) * f
        return fx, fy

def loop(bodies):
    """([Body])

    Never returns; loops through the simulation, updating the
    positions of all the provided bodies.
    """
    timestep = 24*3600  # One day

    step = 1
    while True:
        #commented out update function
        #update_info(step, bodies)
        step += 1

        force = {}
        for body in bodies:
            # Add up all of the forces exerted on 'body'.
            total_fx = total_fy = 0.0
            for other in bodies:
                # Don't calculate the body's attraction to itself lol
                if body is other:
                    continue
                fx, fy = body.attraction_force(other)
                total_fx += fx
                total_fy += fy

            # Record the total force exerted.
            force[body] = (total_fx, total_fy)

        # Update velocities based upon on the force.
        for body in bodies:
            sleep(0.0001)
            fx, fy = force[body]
            body.vx += fx / body.mass * timestep
            body.vy += fy / body.mass * timestep

            # Update positions
            body.px += body.vx * timestep
            body.py += body.vy * timestep
            newpos = vector(body.px, body.py, 0)
            body.model.pos = newpos

def main():
    sun = Body()
    sun.name = 'Sun'
    sun.mass = M_SUN
    sun.model = sphere(pos = vector(0,0,0),
                       radius = R_SUN * SUNSCALE,
                       color = color.yellow)


    mercury = Body()
    mercury.name = 'Mercury'
    mercury.mass = M_MER
    mercury.px = 0.39 * AU
    mercury.vy = -47.87 * 1000
    mercury.model = sphere(pos = vector(mercury.px,0,0),
                           radius = R_MER * SMALLBODYSCALE,
                           color = vec(1,1,1))

    venus = Body()
    venus.name = 'Venus'
    venus.mass = M_VEN
    venus.px = 0.723 * AU
    venus.vy = -35.02 * 1000
    venus.model = sphere(pos = vector(venus.px,0,0),
                         radius = R_VEN * SMALLBODYSCALE,
                         color = color.orange)

    earth = Body()
    earth.name = 'Earth'
    earth.mass = M_EAR
    earth.px = -1 * AU
    earth.vy = 29.783 * 1000
    earth.model = sphere(pos = vector(earth.px,0,0),
                         radius = R_EAR * SMALLBODYSCALE,
                         color = color.blue)

    mars = Body()
    mars.name = 'Mars'
    mars.mass = M_MAR
    mars.px = -1.524 * AU
    mars.vy = 24.1 * 1000
    mars.model = sphere(pos = vector(mars.px,0,0),
                         radius = R_MAR * SMALLBODYSCALE,
                         color = color.red)

    jupiter = Body()
    jupiter.name = 'Jupiter'
    jupiter.mass = M_JUP
    jupiter.px = 5.203 * AU
    jupiter.vy = -13.1 * 1000
    jupiter.model = sphere(pos = vector(jupiter.px,0,0),
                         radius = R_JUP * LARGEBODYSCALE,
                         color = color.green)

    saturn = Body()
    saturn.name = 'Saturn'
    saturn.mass = M_SAT
    saturn.px = 9.539 * AU
    saturn.vy = -9.69 * 1000
    saturn.model = sphere(pos = vector(saturn.px,0,0),
                         radius = R_SAT * LARGEBODYSCALE,
                         color = vec(1,0,1))

    uranus = Body()
    uranus.name = 'Uranus'
    uranus.mass = M_URA
    uranus.px = 19.18 * AU
    uranus.vy = -6.81 * 1000
    uranus.model = sphere(pos = vector(uranus.px,0,0),
                         radius = R_URA * LARGEBODYSCALE,
                         color = color.white)

    # Note, Neptune appears to be quite far away. this makes relative
    # sizes very small
    neptune = Body()
    neptune.name = 'Neptune'
    neptune.mass = M_NEP
    neptune.px = 30.06 * AU
    neptune.vy = -5.43 * 1000
    neptune.model = sphere(pos = vector(neptune.px,0,0),
                         radius = R_NEP * LARGEBODYSCALE,
                         color = color.blue)

    loop([sun, mercury, venus, earth, mars, jupiter, saturn,
          uranus, neptune])
if __name__ == '__main__':
    main()
