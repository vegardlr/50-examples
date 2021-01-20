#!/usr/bin/env python3

import math
import numpy as np
import matplotlib.pyplot as plt
from turtle import *

# The gravitational constant G
G = 6.67428e-11

# Assumed scale: 100 pixels = 1AU.
AU = (149.6e6 * 1000)     # 149.6 million km, in meters.
SCALE = 250 / AU
#timestep = 0.5*3600 # Half hour
timestep = 1*3600 # One hour
#timestep = 12*3600 # Half days
#timestep = 24*3600  # One day
#timestep = 96*3600 # Four days

max_steps = int(365*24*3600/timestep) #One Earth year
max_steps = 1000

class Body(Turtle):
    """Subclass of Turtle representing a gravitationally-acting body.

    Extra attributes:
    mass : mass in kg
    vx, vy: x, y velocities in m/s
    px, py: x, y positions in m
    """
    
    name = 'Body'
    mass = None
    radius = None
    vx = vy = 0.0
    px = py = 0.0
    #Format: step,px,py,vx,vy,fx,fy
    data = None
 #   print("Shape:",data.shape)

    def init(self):
        self.data = np.zeros((max_steps+1,7))


    
    def attraction(self, other):
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
        minimum_distance = self.radius + other.radius

        # Report an error if the distance is zero; otherwise we'll
        # get a ZeroDivisionError exception further down.
        if d <= minimum_distance:
            raise ValueError("Collision between objects %r and %r"
                             % (self.name, other.name))

        # Compute the force of attraction
        f = G * self.mass * other.mass / (d**2)

        # Compute the direction of the force.
        theta = math.atan2(dy, dx)
        fx = math.cos(theta) * f
        fy = math.sin(theta) * f
        return fx, fy

def update_info(step, bodies):
    """(int, [Body])
    
    Displays information about the status of the simulation.
    """
    print('Step #{}'.format(step))
    for body in bodies:
        vabs = math.sqrt(body.vx**2 + body.vy**2)
        s = '{:<8}  Pos.={:>12.8f} {:>12.8f} Vel.={:>12.5f} {:>12.5f} (|{:>12.5f}|)'.format(
            body.name, body.px/AU, body.py/AU, body.vx, body.vy, vabs)
        print(s)
    print()

def init(bodies):
    """([Body])

    Shifts velocities with a half timestep backwards to achieve
    staggered time steps. Reference: 
    Leapfrog method: https://en.wikipedia.org/wiki/Leapfrog_integration
    """
    force = {}
    for body in bodies:
        body.init()
        # Add up all of the forces exerted on 'body'.
        total_fx = total_fy = 0.0
        for other in bodies:
            # Don't calculate the body's attraction to itself
            if body is other:
                continue
            fx, fy = body.attraction(other)
            total_fx += fx
            total_fy += fy

        # Record the total force exerted.
        force[body] = (total_fx, total_fy)

    # Update velocities based upon on the force, half timestep backwards
    for body in bodies:
        fx, fy = force[body]
        body.vx -= fx / body.mass * timestep / 2
        body.vy -= fy / body.mass * timestep / 2
        body.data[0] = np.array([0,body.px,body.py,body.vx,body.vy,fx,fy])
       
def loop(bodies):
    """([Body])

    Never returns; loops through the simulation, updating the
    positions of all the provided bodies.
    """
    update_info(0, bodies)
    init(bodies)


    for body in bodies:
        body.penup()
        body.hideturtle()

    step = 0
    while True:
        #update_info(step, bodies)
        step += 1

        force = {}
        for body in bodies:
            # Add up all of the forces exerted on 'body'.
            total_fx = total_fy = 0.0
            for other in bodies:
                # Don't calculate the body's attraction to itself
                if body is other:
                    continue
                fx, fy = body.attraction(other)
                total_fx += fx
                total_fy += fy

            # Record the total force exerted.
            force[body] = (total_fx, total_fy)
           
        # Update velocities based upon on the force.
        for body in bodies:
            fx, fy = force[body]
            body.vx += fx / body.mass * timestep
            body.vy += fy / body.mass * timestep

            # Update positions
            body.px += body.vx * timestep
            body.py += body.vy * timestep
            body.goto(body.px*SCALE, body.py*SCALE)
            body.dot(3)
            body.data[step] = np.array([step,body.px,body.py,body.vx,body.vy,fx,fy])

           
        if (step >= max_steps): 
            break

def diagnose(bodies):
    """([Body])

    Plots saved data and checks conservation of energy and momentum.
    """
    #data format: step,px,py,vx,vy,fx,fy
    #for body in bodies:
        #print("")
        #print("########################")
        #print(body.name)
        #print("0:---------------")
        #print(body.data[0])
        #print("1:---------------")
        #print(body.data[1])
        #print("2:---------------")
        #print(body.data[2])
        #print("3:---------------")
        #print(body.data[3])
        #print("4:---------------")
        #print(body.data[4])
        #print("5:---------------")
        #print(body.data[5])
        #print("6:---------------")
        #print(body.data[6])
        #print("7:---------------")
        #print(body.data[7])
    #return 0

    
    for body in bodies:
        fig,axs = plt.subplots(3, sharex=True)
        print(body.name)
        step = [d[0] for d in body.data]
        px =   [d[1] for d in body.data]
        py =   [d[2] for d in body.data]
        vx =   [d[3] for d in body.data]
        vy =   [d[4] for d in body.data]
        fx =   [d[5] for d in body.data]
        fy =   [d[6] for d in body.data]
        vabs = [math.sqrt(d[3]**2+d[4]**2) for d in body.data]
        Fabs = [math.sqrt(d[5]**2+d[6]**2) for d in body.data]
        axs[0].plot(vx,'bo')
        axs[0].plot(vy,'r+')
        axs[1].plot(vabs,'g-')
        axs[2].plot(Fabs,'b+')
        fig.suptitle(body.name)
        axs[0].set(ylabel="Vx,Vy")
        axs[1].set(ylabel="|V|")
        axs[2].set(ylabel="|F|")
        plt.show()
    
    Ek_sum = np.zeros(max_steps+1)
    n = len(bodies)
    fig,axs = plt.subplots(n+1, sharex=True)
    i=0
    for body in bodies:
        vx =   [d[3] for d in body.data]
        vy =   [d[4] for d in body.data]
        Ek = [0.5*body.mass*(d[3]**2+d[4]**2) for d in body.data]
        Ek_sum += np.array(Ek)
        axs[i].plot(Ek)
        axs[i].set(ylabel="Ek ("+body.name+")")
        i += 1
    
    axs[n].plot(Ek_sum)
    axs[n].set(ylabel="TOTAL Ek")
    plt.show()


        


def main():
    sun = Body()
    sun.name = 'Sun'
    sun.mass = 1.98892 * 10**30 #kg
    sun.radius = 6.957 * 10**8 #m 
    sun.pencolor('yellow')

    earth = Body()
    earth.name = 'Earth'
    earth.mass = 5.9742 * 10**24
    earth.radius = 6.3781 * 10**6 #m
    earth.px = -1*AU
    earth.vy = -29.783 * 1000            # 29.783 km/sec
    earth.pencolor('blue')

    # Venus parameters taken from
    # http://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
    venus = Body()
    venus.name = 'Venus'
    venus.mass = 4.8685 * 10**24
    venus.radius = 6.051 * 10**6 #m
    venus.px = 0.723 * AU
    venus.vy = 35.02 * 1000
    venus.pencolor('red')

    moon = Body()
    moon.name = 'Moon'
    moon.mass = 7342 * 10**22
    moon.radius = 1.7381 * 10**6
    moon.px = earth.px + 384399000.0 #m
    moon.py = earth.py
    moon.vx = 0.0
    moon.vy = earth.vy + 1022.0 #m/s
    moon.pencolor('gray')

    loop([sun, earth, venus, moon])
    diagnose([sun, earth, venus, moon])


if __name__ == '__main__':
    main()
