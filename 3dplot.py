from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas

points1 = pandas.read_csv('foo1.csv')
points2 = pandas.read_csv('foo2.csv')
points3 = pandas.read_csv('foo3.csv')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x1 = points1['x'].values
y1 = points1['y'].values
z1 = points1['z'].values
x2 = points2['x'].values
y2 = points2['y'].values
z2 = points2['z'].values
x3 = points3['x'].values
y3 = points3['y'].values
z3 = points3['z'].values


ax.plot(x1, y1, z1)
ax.plot(x2, y2, z2)
ax.plot(x3, y3, z3)

ax.scatter([x1[-1],x2[-1],x3[-1]],[y1[-1],y2[-1],y3[-1]],[z1[-1],z2[-1],z3[-1]],"o")
plt.show()