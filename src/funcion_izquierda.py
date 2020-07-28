# export LD_LIBRARY_PATH=/home/user/lab_ws/install/lib:$LD_LIBRARY_PATH

#export PYTHONPATH=/home/user/lab_ws/install/lib/python2.7/site-packages:$PYTHONPATH

import numpy as np
from copy import copy
import rbdl

pi = np.pi

cos=np.cos; sin=np.sin; pi=np.pi


class Robot(object):
    def __init__(self, q0, dq0, ndof, dt):
        self.q = q0    # numpy array (ndof x 1)
        self.dq = dq0  # numpy array (ndof x 1)
        self.M = np.zeros([ndof, ndof])
        self.b = np.zeros(ndof)
        self.dt = dt
        self.robot = rbdl.loadModel('../urdf/robot_gazebo_proyecto_izq.urdf')

    def send_command(self, tau):
        rbdl.CompositeRigidBodyAlgorithm(self.robot, self.q, self.M)
        rbdl.NonlinearEffects(self.robot, self.q, self.dq, self.b)
        ddq = np.linalg.inv(self.M).dot(tau-self.b)
        self.q = self.q + self.dt*self.dq
        self.dq = self.dq + self.dt*ddq

    def read_joint_positions(self):
        return self.q

    def read_joint_velocities(self):
        return self.dq


def DH(d, theta, a, alpha):
    """
    Calcular la matriz de transformacion homogenea asociada con los parametros
    de Denavit-Hartenberg.
    Los valores d, theta, a, alpha son escalares.
    """
    # Escriba aqui la matriz de transformacion homogenea en funcion de los valores de d, theta, a, alpha
    TDH = np.array([[cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta)],
		            [sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
		            [0         ,  sin(alpha)           ,  cos(alpha)           , d           ],
		            [0         ,  0                    ,  0                   , 1           ]])

    T = np.round(TDH,5)
    return T
    
    

def fkine(q):
    """
	
    """
    # Longitudes (en metros)
    l0 = 0.420
    l6 = 0.250
    l7 = 0.06
    l8 = 0.320
    l9 = 0.180
    l10 = 0.15

    Tb =  np.array([[1,0,0,0],[0,1,0,l6],[0,0,1,l0],[0,0,0,1]])
    T01 = DH( l7,    q[0],  l8, 0)
    T12 = DH(  0,    q[1]-np.pi/2,  l9, np.pi)
    T23 = DH( l10+q[2], 0,0,0)



    # Efector final con respecto a la base
    T = Tb.dot(T01).dot(T12).dot(T23)

    return T


def jacobian_position(q, delta=0.0001):
    """
    Jacobiano analitico para la posicion. Retorna una matriz de 3x6 y toma como
    entrada el vector de configuracion articular q=[q1, q2, q3, q4, q5, q6]
    """
    # Crear una matriz 3x6
    J = np.zeros((3,3))
    
    # Transformacion homogenea inicial (usando q)
    T = fkine(q)
    
    # Iteracion para la derivada de cada columna
    for i in xrange(3):
        # Copiar la configuracion articular inicial
        dq = copy(q)

        # Incrementar la articulacion i-esima usando un delta
        dq[i] = dq[i] + delta     

        # Transformacion homogenea luego del incremento (q+delta)
        Td = fkine(dq) 

        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        J[0,i] = (1/delta)*(Td[0,3]-T[0,3])
        J[1,i] = (1/delta)*(Td[1,3]-T[1,3])
        J[2,i] = (1/delta)*(Td[2,3]-T[2,3])

    J = np.round(J,3)
    return J





def rot2quat(R):
    """
    Convertir una matriz de rotacion en un cuaternion

    Entrada:
      R -- Matriz de rotacion
    Salida:
      Q -- Cuaternion [ew, ex, ey, ez]

    """
    dEpsilon = 1e-6;
    quat = 4*[0.,]

    quat[0] = 0.5*np.sqrt(R[0,0]+R[1,1]+R[2,2]+1.0)
    if ( np.fabs(R[0,0]-R[1,1]-R[2,2]+1.0) < dEpsilon ):
        quat[1] = 0.0
    else:
        quat[1] = 0.5*np.sign(R[2,1]-R[1,2])*np.sqrt(R[0,0]-R[1,1]-R[2,2]+1.0)
    if ( np.fabs(R[1,1]-R[2,2]-R[0,0]+1.0) < dEpsilon ):
        quat[2] = 0.0
    else:
        quat[2] = 0.5*np.sign(R[0,2]-R[2,0])*np.sqrt(R[1,1]-R[2,2]-R[0,0]+1.0)
    if ( np.fabs(R[2,2]-R[0,0]-R[1,1]+1.0) < dEpsilon ):
        quat[3] = 0.0
    else:
        quat[3] = 0.5*np.sign(R[1,0]-R[0,1])*np.sqrt(R[2,2]-R[0,0]-R[1,1]+1.0)

    return np.array(quat)


def TF2xyzquat(T):
    """
    Convert a homogeneous transformation matrix into the a vector containing the
    pose of the robot.

    Input:
      T -- A homogeneous transformation
    Output:
      X -- A pose vector in the format [x y z ew ex ey ez], donde la first part
           is Cartesian coordinates and the last part is a quaternion
    """
    quat = rot2quat(T[0:3,0:3])
    res = [T[0,3], T[1,3], T[2,3], quat[0], quat[1], quat[2], quat[3]]
    return np.array(res)



def skew(w):
    R = np.zeros([3,3])
    R[0,1] = -w[2]; R[0,2] = w[1];
    R[1,0] = w[2];  R[1,2] = -w[0];
    R[2,0] = -w[1]; R[2,1] = w[0];
    return R


def ikine(xdes, q0):
    """

    """
    epsilon  = 0.001
    max_iter = 100000
    delta    = 0.00001

    q  = copy(q0)
    for i in range(max_iter):
        # Jacobiano
        J = jacobian_position(q, delta=0.0001)
        # Cinematica directa
        F = fkine(q)
        f = F[0:3,3]
        e = xdes-f
        q = q + np.dot(np.linalg.pinv(J), e)
        # Condicion de termino
        if (np.linalg.norm(e) < epsilon):
            break       
    return q
