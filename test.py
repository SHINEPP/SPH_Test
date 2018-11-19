# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from SPHSystem import SPHSystem

sph = SPHSystem()
sph.set_running_flag(1)
count = sph.setup_system()

x_data = [0] * count
y_data = [0] * count
z_data = [0] * count


def update_inner():
    sph.update_particle()
    mem = sph.get_mem()
    for i in range(0, count):
        p = mem[i]
        x_data[i] = p.position[0]
        y_data[i] = p.position[1]
        z_data[i] = p.position[2]


def update_graph(num):
    sph.update_particle()
    mem = sph.get_mem()

    for i in range(0, count):
        p = mem[i]
        x_data[i] = p.position[0]
        y_data[i] = p.position[2]
        z_data[i] = p.position[1]

    graph.set_data(x_data, y_data)
    graph.set_3d_properties(z_data)
    title.set_text('SPH Test, time={}'.format(num))
    return title, graph,


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('SPH Test')
ax.set_xlabel('x')
ax.set_ylabel('y')

graph, = ax.plot([], [], [], linestyle="", marker=".")

ani = matplotlib.animation.FuncAnimation(fig, update_graph, 10, interval=100, blit=True)

plt.show()
