import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import pandas

points1 = pandas.read_csv('foo1.csv')
points2 = pandas.read_csv('foo2.csv')
points3 = pandas.read_csv('foo3.csv')
points4 = pandas.read_csv('foo4.csv')
points5 = pandas.read_csv('foo5.csv')

x1 = points1['x'].values
y1 = points1['y'].values
z1 = points1['z'].values
x2 = points2['x'].values
y2 = points2['y'].values
z2 = points2['z'].values
x3 = points3['x'].values
y3 = points3['y'].values
z3 = points3['z'].values
x4 = points4['x'].values
y4 = points4['y'].values
z4 = points4['z'].values
x5 = points5['x'].values
y5 = points5['y'].values
z5 = points5['z'].values

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1, y1, z1)
ax.plot(x2, y2, z2)
ax.plot(x3, y3, z3)
ax.plot(x4, y4, z4)
ax.plot(x5, y5, z5)
ax.scatter([x1[-1],x2[-1],x3[-1],x4[-1],x5[-1]],[y1[-1],y2[-1],y3[-1],y4[-1],y5[-1]],[z1[-1],z2[-1],z3[-1],z4[-1],z5[-1]],"o")

def anim_func(i):
    ax.plot(x1[:50*i:5], y1[:50*i:5], z1[:50*i:5])
    ax.plot(x2[:50*i:5], y2[:50*i:5], z2[:50*i:5])
    ax.plot(x3[:50*i:5], y3[:50*i:5], z3[:50*i:5])
    ax.plot(x4[:50*i:5], y4[:50*i:5], z4[:50*i:5])
    ax.plot(x5[:50*i:5], y5[:50*i:5], z5[:50*i:5])
    #ax.scatter([x1[50*i],x2[50*i],x3[50*i],x4[-1],x5[-1]],[y1[-1],y2[-1],y3[-1],y4[-1],y5[-1]],[z1[-1],z2[-1],z3[-1],z4[-1],z5[-1]],"o")
animation = FuncAnimation(fig, anim_func, frames =  200,interval = 20)
animation.save("trial2.mp4", dpi=600)