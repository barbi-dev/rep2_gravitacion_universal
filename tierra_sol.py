from vpython import sphere, vector, rate, color, scene, textures
import numpy as np

G=4*np.pi**2   # ua
softening=1e-9     
dt=0.0005       #h: tamaño de paso 
Tmax=2           
frame_rate=200     

m_sun=1
m_earth=3.003e-6
masses=np.array([m_sun, m_earth])
# orbitales (elipse)
a = 1        
e = 0.0167    
r0 = a*(1 - e)           
v0 = np.sqrt(G*m_sun*(2/r0-1/a))
# sol [0], tierra [1]
pos = np.zeros((2, 3), dtype=float)
vel = np.zeros((2, 3), dtype=float)
pos[0] = np.array([0, 0, 0])     # sol
pos[1] = np.array([0, r0, 0])      # tierra

vel_earth=np.array([-v0, 0, 0])
vel_sun=-(m_earth/m_sun)*vel_earth   
vel[0]=vel_sun
vel[1]=vel_earth

# VPython
scene.background=color.black
scene.height=540
scene.width=860
scene.forward=vector(-1,-0.2,-0.3)
scene.range=0.9
sun_sphere=sphere(pos=vector(*pos[0]), radius=0.2, color=color.yellow, emissive=True, make_trail=False, texture=textures.stucco)
earth_sphere=sphere(pos=vector(*pos[1]), radius=0.08, make_trail=True, retain=300, texture=textures.earth)

# aceleraciones 
def accelerations(positions, masses):
    N=positions.shape[0]
    acc=np.zeros_like(positions)
    for i in range(N):
        ai=np.zeros(3)
        for j in range(N):
            if i==j:
                continue
            rvec=positions[j]-positions[i]
            dist=np.linalg.norm(rvec)
            denom=(dist**3)+softening
            ai+= G*masses[j]*rvec/denom
        acc[i]=ai
    return acc
# RK4 
def rk4(pos, vel, masses, dt):
    # k1
    a1=accelerations(pos, masses)
    k1r=vel
    k1v=a1
    # k2
    pos_k2=pos+0.5*dt*k1r
    vel_k2=vel+0.5*dt*k1v
    a2=accelerations(pos_k2, masses)
    k2r=vel_k2
    k2v=a2
    # k3
    pos_k3=pos+0.5*dt*k2r
    vel_k3=vel+0.5*dt*k2v
    a3=accelerations(pos_k3, masses)
    k3r=vel_k3
    k3v=a3
    # k4
    pos_k4=pos+dt*k3r
    vel_k4=vel+dt*k3v
    a4=accelerations(pos_k4, masses)
    k4r=vel_k4
    k4v=a4

    pos_new=pos+(dt/6)*(k1r+2*k2r+2*k3r+k4r)
    vel_new=vel+(dt/6)*(k1v+2*k2v+2*k3v+k4v)

    return pos_new, vel_new

# simulación
t=0
while t < Tmax:
    rate(frame_rate)
    pos,vel=rk4(pos, vel, masses, dt)

    sun_sphere.pos=vector(*pos[0])
    earth_sphere.pos=vector(*pos[1])
    earth_sphere.rotate(angle=0.01,axis=vector(0,1,0))
    scene.camera.follow(earth_sphere)
    t += dt
#IMPORTANTE
#*pos[0] convierte array([x0, y0, z0]) en tres argumentos separados x0, y0, z0
#para poder usarlos con el comando vector()