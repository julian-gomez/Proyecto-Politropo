import numpy as np
import matplotlib.pyplot as plt

datos = np.genfromtxt("casoAnisotropo.txt", delimiter = " ", usecols = (0,1))
plt.figure()
plt.plot(datos[:,0],datos[:,1])
plt.show()
