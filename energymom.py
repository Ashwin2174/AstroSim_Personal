import numpy as np
import matplotlib.pyplot  as plt

e,px,py,pz = np.loadtxt("foo1.csv", unpack=True, delimiter=",")

plt.title("Energy of the System")
plt.plot([i for i in range(len(e))],e)
plt.show()
