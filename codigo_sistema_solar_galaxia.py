from vpython import vector, sphere, color, rate, scene, curve
import math
import random

TWOPI = 2*math.pi
#PLANTEO ECUACION DE KEPLER    M= E-e*sin(E)
#M anomalía media: tiempo convertido en ángulo, cuanto avanza planeta con omega constante
#M es angulo ficticio, como si el planeta se moviera en un círculo.
#E anomalía excéntrica: angulo aux. permite pasar de la órbita circular a elíptica.
def kepler_E(M, e, tol=1e-9, maxiter=50):
    # normalizar M para que quede entre 0 y 2pi
    M = M % TWOPI
    #valor inicial de E para newton raphson
    #si orbita no exentrica E=M coonverge rapido, si no se lanza guess=pi
    if e < 0.8:
        E = M
    else:
        E = math.pi
    for _ in range(maxiter):
        f = E - e * math.sin(E) - M
        df = 1 - e * math.cos(E)
        dE = -f / df
        E += dE
        if abs(dE) < tol:
            break
        #devuelve la anomalía excéntrica E para calcular la anomalía verdadera v (pos angular real en orbita eliptica)
    return E

#ajustamos angulo del planeta en su orbita tal que el perihelio se ubique en angulo (omega)
def rotate_about_z(v, ang):
    c=math.cos(ang);s= math.sin(ang)
    #uso formula para rotacion alrededor de z: x'=xcostheta-ysintheta , y'=xsintheta+ycostheta, z'=z
    return vector(c*v.x - s*v.y, s*v.x + c*v.y, v.z)

#rota el plano orbital para darle inclinacion respecto al plano de referencia (se le llama nodo ascendente)
def rotate_about_x(v, ang):
    c = math.cos(ang); s = math.sin(ang)
    return vector(v.x, c*v.y - s*v.z, s*v.y + c*v.z)

# devuelve la posición exacta del planeta en el sistema inercial
#Rota toda la órbita para ubicarla en el espacio real
def orbital_to_inertial(r_orb, Omega, inclination, omega):
    v = rotate_about_z(r_orb, omega)    
    v = rotate_about_x(v, inclination)        
    v = rotate_about_z(v, Omega)    
    return v

#MOVIMIENTO GALACTICO
#POS DEL SOL EN FUNCION DEL TIEMPO
#amp y freq muestran pequeña oscilacion tipo onda; base es la pos inicial del sol; v velocidad
def sun_galactic_position(t, params):
    v = params['v']        
    base = params['sun0']  
    amp = params['amp']    
    freq = params['freq'] ###ciclos por año
    #movimiento rectilíneo uniforme del sol
    pos = base + vector(v[0]*t, v[1]*t, v[2]*t)
    # añade pequeño movimiento oscilatorio
    pos += vector(amp * math.cos(TWOPI*freq*t), amp * math.sin(TWOPI*freq*t), amp * math.sin(0.5*TWOPI*freq*t))
    return pos

# Parámetros
M_central = 1 # masa del Sol en masas solares
scale_visual = 1 # escala multiplicativa para posiciones

# Planetas: lista con parámetros
# a (AU), e, incl (rad), Omega (rad), omega (rad), phase (rad initial M0), color, size (visual radius)
planets = [
    {"name":"Mercurio", "a":0.39, "e":0.205, "inc":0.0, "Omega":0.0, "omega":0.0, "M0":random.random()*TWOPI, "col":color.gray(0.6), "r":0.04},
    {"name":"Venus",   "a":0.72, "e":0.0067, "inc":0.0, "Omega":0.0, "omega":0.3, "M0":random.random()*TWOPI, "col":color.orange, "r":0.06},
    {"name":"Tierra",  "a":1.00, "e":0.0167, "inc":0.0, "Omega":0.0, "omega":0.5, "M0":0.0, "col":color.cyan, "r":0.08},
    {"name":"Marte",   "a":1.25, "e":0.0934, "inc":0.0, "Omega":0.0, "omega":1.0, "M0":random.random()*TWOPI, "col":color.red, "r":0.04},
    
]

# Calcular períodos y velocidades angulares (unidades: years)
for p in planets:
    a = p["a"]
    # Para unidades AU & years & masas solares: periodo^2 = a^3 / M_central
    period = math.sqrt(a**3/M_central)
    p["period"] = period
    p["n"] = TWOPI / period

# Crear escena VPython
scene.title = "Sistema planetario + movimiento galáctico"
scene.width = 860
scene.height = 540
scene.background = color.black
scene.range = 1.5  
# Sol (centrado dinámicamente por la función sun_galactic_position)
sun = sphere(pos=vector(0,0,0), radius=0.2, color=color.yellow, emissive=True, make_trail=True, retain=150)
sun_label_offset = vector(0.5, 0.5, 0)

# lista de objetos para planetas y rastro
for p in planets:
    p["body"] = sphere(pos=vector(0,0,0), radius=p["r"], color=p["col"], make_trail=True, retain=100)
    p["trail_curve"] = curve(color=p["col"])

# Parámetros del movimiento galáctico del Sol
galaxy_params = {
    "v": (6.0, 0.5, 0.0),
    "sun0": vector(0, 0.0, 0.0),  
    "amp": 0.9,  
    "freq": 0.1 
}

# Ajustes temporales
t = 0.0
dt = 0.005   # years por paso (aprox 0.73 días); ajusta para velocidad/precisión
max_steps = 2000

# Para ver órbitas, también podemos dibujar la órbita geométrica
def draw_orbit_path(p, steps=360):
    pts = []
    a = p["a"]; e = p["e"]
    b = a * math.sqrt(1 - e*e)
    for i in range(steps):
        theta = TWOPI * i / steps
        # param enorb: usar excentric anomaly para mu uniform: elipse geométrica
        x = a * (math.cos(theta) - e)
        y = b * math.sin(theta)
        ro = orbital_to_inertial(vector(x, y, 0), p["Omega"], p["inc"], p["omega"])
        pts.append(ro * scale_visual + sun_galactic_position(0.0, galaxy_params))  # se centra inicialmente
    return curve(points=pts, radius=0.002, color=p["col"], opacity=0.15)

# Dibuja las elipses iniciales
orbit_curves = []
for p in planets:
    orbit_curves.append(draw_orbit_path(p))

# ---------- MAIN LOOP ----------

while t < max_steps * dt:
    rate(200)  
    sun_pos = sun_galactic_position(t, galaxy_params)
    sun.pos = sun_pos * scale_visual

    # actualiza cada planeta
    for p in planets:
        # anomalia media
        M = p["M0"] + p["n"] * t
        # resolver Kepler
        E = kepler_E(M, p["e"])
        # posición en el plano orbital (parametrico)
        a = p["a"]
        x_orb = a * (math.cos(E) - p["e"])
        y_orb = a * math.sqrt(max(0.0, 1.0 - p["e"]**2)) * math.sin(E)
        r_orb = vector(x_orb, y_orb, 0.0)
        # convertir a coordenadas inerciales (rotaciones Omega, i, omega)
        r_inertial = orbital_to_inertial(r_orb, p["Omega"], p["inc"], p["omega"])
        # pos final = pos del Sol + r_inertial
        pos = sun_pos * scale_visual + r_inertial * scale_visual
        scene.center = sun.pos                  
        scene.forward = vector(-1, -0.5, -0.5)     
        scene.up = vector(0, 0, 1) 
        # actualizar cuerpo
        p["body"].pos = pos
        # actualización de trail
        p["trail_curve"].append(pos)

    t += dt

