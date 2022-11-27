import numpy as np
import matplotlib.pylab as plt

datos0 = np.loadtxt("datos_runge_kutta0.txt", delimiter= ",")
datos1 = np.loadtxt("datos_runge_kutta1.txt", delimiter= ",")
datos = np.loadtxt("datos_runge_kutta.txt", delimiter= ",")

g_0 = plt.figure()
circle1 = plt.Circle((0, 0), 6371, color='r')
plt.gca().add_patch(circle1)
plt.plot(datos0[:,1]*np.cos(datos0[:,2]),datos0[:,1]*np.sin(datos0[:,2]))
plt.plot(datos1[:,1]*np.cos(datos1[:,2]),datos1[:,1]*np.sin(datos1[:,2]))
###plt.plot(datos[:,1]*np.cos(datos[:,2]),datos[:,1]*np.sin(datos[:,2]))
###plt.plot(datos[:,0],datos[:,1]*np.sin(datos[:,2]))
###plt.xlim([6371,6390])
###plt.plot(datos[:,0],datos[:,4])
###plt.plot(datos[:,0], datos[:,1])
plt.savefig("Resorte.png")