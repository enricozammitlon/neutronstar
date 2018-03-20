import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon

fig, ax = plt.subplots()
ax.add_patch(Wedge((0, 0), 12, 0, 90))
#Use adjustable='box-forced' to make the plot area square-shaped as well.
#ax.set_aspect('equal', adjustable='datalim')
ax.plot()   #Causes an autoscale update.
plt.show()
