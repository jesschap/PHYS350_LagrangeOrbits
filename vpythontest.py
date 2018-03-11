# pylama:ignore=W391,E231,E251,E303,E221,E272

import math
from vpython import *

# Gravitational constant
G = 6.67428e-11

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


#class Body():
#    """
#    mass : mass in kg
#    vx, vy: x, y velocities in m/s
#    px, py: x, y positions in m
#    """
#
#    name = None
#    mass = None
#    vx = vy = 0.0
#    px = py = 0.0
#    model = None
#
#
#    # This function computes the force between one body and the other, and then
#    # returns them in the order fx, fy
#    def attraction_force(self, other):
#        """(Body): (fx, fy)
#
#        Returns the force exerted upon this body by the other body.
#        """
#        # Report an error if the other object is the same as this one.
#        print("ehehehhe")
#        if self is other:
#            raise ValueError("Attraction of object %r to itself requested"
#                             % self.name)
#
#        # Compute the distance of the other body.
#        sx, sy = self.px, self.py
#        ox, oy = other.px, other.py
#        dx = (ox-sx)
#        dy = (oy-sy)
#        d = math.sqrt(dx**2 + dy**2)
#
#        # Report an error if the distance is zero cuz u
#        # get a ZeroDivisionError exception further down.
#        if d == 0:
#            raise ValueError("Collision between objects %r and %r"
#                             % (self.name, other.name))
#
#        # Compute the force of attraction
#        f = G * self.mass * other.mass / (d**2)
#
#        # Compute the direction of the force.
#        theta = math.atan2(dy, dx)
#        fx = math.cos(theta) * f
#        fy = math.sin(theta) * f
#        return fx, fy
#
#def loop(bodies):
#    """([Body])
#
#    Never returns; loops through the simulation, updating the
#    positions of all the provided bodies.
#    """
#    timestep = 24*3600  # One day
#    print("Ijsiddjijds'm here")
#
#    for body in bodies:
#        print("Hello I am " + body.name)
#
#    step = 1
#    while True:
#        #commented out update function
#        #update_info(step, bodies)
#        step += 1
#
#        force = {}
#        for body in bodies:
#            # Add up all of the forces exerted on 'body'.
#            total_fx = total_fy = 0.0
#            for other in bodies:
#                # Don't calculate the body's attraction to itself lol
#                if body is other:
#                    continue
#                fx, fy = body.attraction(other)
#                total_fx += fx
#                total_fy += fy
#
#            # Record the total force exerted.
#            force[body] = (total_fx, total_fy)
#
#        # Update velocities based upon on the force.
#        for body in bodies:
#            fx, fy = force[body]
#            body.vx += fx / body.mass * timestep
#            body.vy += fy / body.mass * timestep
#
#            # Update positions
#            body.px += body.vx * timestep
#            body.py += body.vy * timestep
#            newpos = vector(body.px, body.py, 0)
#            body.model.pos = newpos

def main():
    print("Hello my friend")
#    sun = Body()
#    sun.name = 'Sun'
#    sun.mass = 1.98892 * 10**30
#    sun.px = 0
#    sun.vy = 0
#    sun.model = sphere(pos = vector(0,0,0), radius = 0.5, color = color.yellow)
#
#    print("Hello my friendssss")
#    earth = Body()
#    earth.name = 'Earth'
#    earth.mass = 5.9742 * 10**24
#    earth.px = -1*AU
#    earth.vy = 29.783 * 1000            # 29.783 km/sec
#    earth.model = sphere(pos = vector(10,0,0), radius = 0.5, color = color.blue)
#
#    loop([sun, earth, venus, mars, jupiter, ast])

