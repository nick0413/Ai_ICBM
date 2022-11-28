import numpy as np
import matplotlib.pylab as plt


def graficar(num):
  g_0 = plt.figure()


  circle1 = plt.Circle((0, 0), 6371, color='b')
  #plt.gca().add_patch(circle1)
  for i in range(50):
    d=np.loadtxt(f"datos_runge_kutta{i}.txt",delimiter= ",")
    plt.plot(d[:,1]*np.cos(d[:,2]),d[:,1]*np.sin(d[:,2]))

  plt.savefig(f"trayectorias\grafica{num}.png")
  plt.close()