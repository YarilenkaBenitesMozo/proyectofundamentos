import rbdl
import numpy as np


# Lectura del modelo del robot a partir de URDF (parsing)
modelo = rbdl.loadModel('../urdf/robot_gazebo_proyecto.urdf')
# Grados de libertad
ndof = modelo.q_size

print(ndof)
# Configuracion articular
q = np.array([0.0, -1.0, 1.7, -2.2, -1.6, 0.0])
# Velocidad articular
dq = np.array([0., 0., 0., 0., 0., 0.])
# Aceleracion articular
ddq = np.array([0., 0., 0., 0., 0., 0.])

# Arrays numpy
zeros = np.zeros(ndof)          # Vector de ceros
tau   = np.zeros(ndof)          # Para torque
g     = np.zeros(ndof)          # Para la gravedad
c     = np.zeros(ndof)          # Para el vector de Coriolis+centrifuga
M     = np.zeros([ndof, ndof])  # Para la matriz de inercia
e     = np.eye(3)               # Vector identidad

# Torque dada la configuracion del robot
rbdl.InverseDynamics(modelo, q, dq, ddq, tau)
rbdl.CompositeRigidBodyAlgorithm(modelo, q, M)
print(M)
print(tau)
