import numpy as np
import matplotlib.pyplot as plt

datos = np.genfromtxt("masas.txt", delimiter = " ", usecols = (0,1))

plt.figure()
plt.scatter(datos[:,0],datos[:,1],label = "m = 3.0, alpha = beta = 1.0")
plt.title("Indice M_iso/M_ani variando n")
plt.legend()
plt.savefig("Indice-M_iso-M_ani-variando-n.jpg")
