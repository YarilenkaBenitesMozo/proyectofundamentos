#!/usr/bin/env python

import rospy
from sensor_msgs.msg import JointState
from markers import *
from funcion_izquierda import *
from roslib import packages

import rbdl
import numpy as np

rospy.init_node("control_pdg")
pub = rospy.Publisher('joint_states', JointState, queue_size=1000)
bmarker_actual  = BallMarker(color['RED'])
bmarker_deseado = BallMarker(color['GREEN'])
# Archivos donde se almacenara los datos
fqact = open("/tmp/qactual.dat", "w")
fqdes = open("/tmp/qdeseado.dat", "w")
fxact = open("/tmp/xactual.dat", "w")
fxdes = open("/tmp/xdeseado.dat", "w")

# Nombres de las articulaciones
jnames = ['joint_1I', 'joint_2I','joint_3I']
# Objeto (mensaje) de tipo JointState
jstate = JointState()
# Valores del mensaje
jstate.header.stamp = rospy.Time.now()
jstate.name = jnames

# =============================================================
# Configuracion articular inicial (en radianes)
q = np.array([0., 0., 0.])
# Velocidad inicial
dq = np.array([0., 0., 0.])
# Aceleracion inicial
ddq = np.array([0., 0., 0.])
# Configuracion articular deseada
qdes = np.array([0.85,0.0,0.06])
# Velocidad articular deseada
dqdes = np.array([0.0, 0.0, 0.0])
# Aceleracion articular deseada
ddqdes = np.array([0.0, 0.0, 0.0])
# =============================================================

# Posicion resultante de la configuracion articular deseada
T = fkine(qdes)
xdes = T[0:3,3]
# Copiar la configuracion articular en el mensaje a ser publicado
jstate.position = q
pub.publish(jstate)

# Modelo RBDL
modelo = rbdl.loadModel('../urdf/robot_gazebo_proyecto_izq.urdf')
ndof   = modelo.q_size     # Grados de libertad
zeros = np.zeros(ndof)     # Vector de ceros

# Frecuencia del envio (en Hz)
freq = 250
dt = 1.0/freq
rate = rospy.Rate(freq)

# Simulador dinamico del robot
robot = Robot(q, dq, ndof, dt)

# Bucle de ejecucion continua
t = 0.0
u = np.zeros(ndof)   # Reemplazar por la ley de control
M = np.zeros([ndof, ndof])
g = np.zeros(ndof)
zeros = np.zeros(ndof)
c = np.zeros(ndof)



# Se definen las ganancias del controlador
valores = 5.0*np.array([1.0, 1.0, 1.0])
Kp = np.diag(valores)
Kd = 2*np.sqrt(Kp)

while not rospy.is_shutdown():

    # Leer valores del simulador
    q  = robot.read_joint_positions()
    dq = robot.read_joint_velocities()
    # Posicion actual del efector final
    x = fkine(q)[0:3,3]
    # Tiempo actual (necesario como indicador para ROS)
    jstate.header.stamp = rospy.Time.now()

    # Almacenamiento de datos


    # ----------------------------
    # Control dinamico (COMPLETAR)
    # ----------------------------

    #VECTOR MASA INERCIAL
    rbdl.CompositeRigidBodyAlgorithm(modelo, q, M)

    #VECTOR GRAVEDAD
    rbdl.NonlinearEffects(modelo, q, zeros, g)

    #VECTOR CORIOLISIS
    rbdl.InverseDynamics(modelo, q, dq, zeros, c)
    c = c-g

    #VECTOR f, CONCATENAMOS VECTOR q Y dq
    f = np.hstack([q,dq])
    
    #VALOR DE y
    y = (ddqdes - ddq )  + Kd.dot(dqdes - f[3:6]) +  Kp.dot(qdes - f[0:3])


    #LEY DE CONTROL
    u = M.dot(y) + c +g  


    #SE ASIGNA NUEVO VECTOR ddq 
    ddq = np.linalg.inv(M).dot(u -g - c )## ace


    #SE OPTIENE dq anterior
    dq = f[3:6]

    #SE OPTIENE vector de ESPACIO DE ESTADOS
    dX = np.hstack([dq,ddq])

    #INTEGRACION NUMERICA
    f = f + dt*dX

    #SE ASIGNA NUEVO VECTOR q
    q = f[0:3]               ### pos

    #SE ASIGNA NUEVO VECTOR dq
    dq = f[3:6]             ### vel



    #Tiempo
    t = t+dt
    #VECTOR FUNCION DE ERROR
    e = qdes - f[0:3]
    print q
    if np.linalg.norm(e)<0.001:
        break

    # Simulacion del robot
    robot.send_command(u)

    # Publicacion del mensaje
    jstate.position = q
    pub.publish(jstate)
    bmarker_deseado.xyz(xdes)
    bmarker_actual.xyz(x)
    t = t+dt
    # Esperar hasta la siguiente  iteracion
    rate.sleep()

fqact.close()
fqdes.close()
fxact.close()
fxdes.close()
