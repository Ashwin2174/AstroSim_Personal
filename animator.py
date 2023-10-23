import matplotlib.pyplot as plt
import csv
from matplotlib.animation import FuncAnimation

x1,y1=[],[]
x2,y2=[],[]
x3,y3=[],[]
with open('trial1.csv', mode ='r')as file:
  csvFile = csv.reader(file)
  for lines in csvFile:
    x1.append(float(lines[0]))
    y1.append(float(lines[1]))
    x2.append(float(lines[2]))
    y2.append(float(lines[3]))
    x3.append(float(lines[4]))
    y3.append(float(lines[5]))
plt.style.use('dark_background')
fig = plt.figure()
ax = plt.axes()
ax.set_ylim(-15,25)
ax.set_xlim(-15,25)
ax.set_aspect(1)
p1 = plt.Circle((0,0), 1, color="blue")
ax.add_artist(p1)
p2 = plt.Circle((0,8), 0.1, color="orange")
ax.add_artist(p2)
p3 = plt.Circle((10,0), 0.1, color="orange")
ax.add_artist(p3)
# plt.plot(x1,y1)
# plt.plot(x2,y2)
# plt.plot(x3,y3)
# plt.show()
plt.rcParams["figure.dpi"] = 400
def anim_func(i):
  if i>=1000:
    ax.plot(x1[:i-1000], y1[:i-1000],color="blue")
    ax.plot(x2[:i-1000], y2[:i-1000],color="orange")
    ax.plot(x3[:i-1000], y3[:i-1000],color="orange")
    p1.center = x1[i-1000],y1[i-1000]
    p2.center = x2[i-1000],y2[i-1000]
    p3.center = x3[i-1000],y3[i-1000]
anim = FuncAnimation(fig, anim_func, frames =  2100,interval = 0.2)
plt.show()

