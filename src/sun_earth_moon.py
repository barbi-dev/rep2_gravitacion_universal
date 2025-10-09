from vpython import sphere, vector, rate, color, textures, curve, scene
import numpy as np

# datos
G=4*np.pi**2
M_sun=1
M_earth=3.003e-6
M_moon=3.69e-8
masses=np.array([M_sun, M_earth,M_moon])
dt=0.0005
steps =1500
#condiciones iniciales
def initial_conditions():
    # posición del CM Tierra-Luna
    r_cm_TL= np.array([1,0])   
    # separación Tierra-Luna en UA
    d_moon=0.00257  
    # posiciones relativas al CM Tierra-Luna 
    r_earth_rel= -d_moon*M_moon/(M_earth + M_moon)
    r_moon_rel=d_moon*M_earth/(M_earth + M_moon)

    r_sun = np.array([0, 0])
    r_earth =r_cm_TL + np.array([r_earth_rel,0])
    r_moon  =r_cm_TL + np.array([r_moon_rel, 0])

    # Velocidad del CM Tierra-Luna alrededor del Sol
    v_cm=np.array([0,2*np.pi]) 
    # velocidad relativa Tierra-luna alrededor del CM 
    v_rel = np.sqrt(G*(M_earth + M_moon)/d_moon)

    # repartimos las velocidades relativas: ambos se mueven en la misma dirección (y positiva)
    v_earth=v_cm + np.array([0,-v_rel*(M_moon/(M_earth + M_moon))])
    v_moon =v_cm + np.array([0,v_rel*(M_earth/(M_earth + M_moon))])

    #Sol con velocidad tal que momento total sea ~0
    v_sun = -(M_earth*v_earth + M_moon*v_moon)/M_sun
    #vstack apila verticalmente los arrays de entrada: todas las filas en una sola matriz
    positions=np.vstack([r_sun,r_earth,r_moon])
    velocities=np.vstack([v_sun,v_earth,v_moon])
    

    return positions, velocities

positions, velocities = initial_conditions()

def compute_accelerations(positions,masses,softening=1e-4):

    # .shape obtiene el número de filas de un objeto array
    N=positions.shape[0]
    #devuelve un nuevo array NumPy que tiene la misma forma y 
    #tipo de datos que un array de entrada dado, con ceros
    acc = np.zeros_like(positions)

    for i in range(N):
        #rij contiene vectores de distancia entre pos de referencia y las demás pos
        rij = positions - positions[i]   # vectores r_j - r_i
        #axis=1): Suma los elementos a lo largo del eje 1
        #el eje de las columnas, es decir, las coordenadas x, y, z.
        dist2 = np.sum(rij**2, axis=1) + softening**2
        inv_dist3 = dist2**-1.5
        inv_dist3[i] = 0.0  # evita división por 0 con sí mismo
        #toma la transpuesta para poder operar las matrices
        acc[i] = G*np.sum((rij.T*(masses * inv_dist3)).T, axis=0)
    return acc

def rk4(positions, velocities, masses, dt):
    def derivative(r,v):
        return v, compute_accelerations(r,masses)

    r1,v1 = positions,velocities
    k1_r,k1_v = derivative(r1,v1)

    r2 = r1 + 0.5*dt*k1_r
    v2 = v1 + 0.5*dt*k1_v
    k2_r, k2_v = derivative(r2,v2)

    r3 = r1 + 0.5*dt*k2_r
    v3 = v1 + 0.5*dt*k2_v
    k3_r, k3_v = derivative(r3,v3)

    r4 = r1 + dt*k3_r
    v4 = v1 + dt*k3_v
    k4_r, k4_v = derivative(r4,v4)

    r_new = r1 + (dt/6)*(k1_r+2*k2_r+2*k3_r+k4_r)
    v_new = v1 + (dt/6)*(k1_v+2*k2_v+2*k3_v+k4_v)

    return r_new, v_new

#cuerpos para la escena
sun = sphere(pos=vector(0,0,0), radius=0.06, color=color.yellow, emissive=True, texture=textures.stucco)
earth = sphere(pos=vector(positions[1][0], positions[1][1], 0),radius=0.0006, texture=textures.earth, make_trail=True, retain=2000)
moon = sphere(pos=vector(positions[2][0], positions[2][1], 0),radius=0.0004, color=color.white, retain=500)
traj = curve(color=color.white)

# simulacion
for i in range(steps):
    rate(60)  
    positions, velocities = rk4(positions, velocities,masses, dt)
    
    # Actualizar posición Tierra
    earth.pos = vector(positions[1][0], positions[1][1],0)
    traj.append(pos=earth.pos)
    moon.pos = vector(positions[2][0], positions[2][1],0)
    # Rotación sobre su eje
    earth.rotate(angle=0.1, axis=vector(0,0,1))
    scene.camera.follow(moon)
