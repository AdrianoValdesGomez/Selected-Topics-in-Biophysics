

# Esta funcion genera, dado un numero, una cadena de una longitud fija (4)
# con ceros a la izquierda. Sirve para poder hacer animaciones usando imagenes
# enumeradas.

def nombre(s):
    diferencia = 4 - len(str(s))
    ceros = '' 
    for i in range(diferencia):
        ceros = ceros + '0'
    variable = ceros + str(s)
    return variable

# Funcion genera una distribucion inicial X0 en una malla
# de dimensiones (d + 2) x (d  + 2)

def dist_uni(d,n):
    X = np.zeros((d + 2,d + 2))
    for i in range(n):
        i,j = np.random.randint(1,d + 1), np.random.randint(1,d + 1)
        X[Si][i][j] = X[Si][i][j] + 1
    return X

# Funcion genera una distribucion circular inicial X0 en una malla
# de dimensiones (d + 2) x (d  + 2)

def dist_circ(centro, radio, n, d):
    X0 = np.zeros((d + 2, d + 2))
    for i in range(1,d + 1):
        for j in range(1,d + 1):
            rs = (i - centro[0])**2 + (j - centro[1])**2
            if rs < radio**2:
                X0[i,j] = n
    return X0



# Funcion que actualiza las propensiones usando el vector de concentracion
# X. Es parte medular del algoritmo de Gillespie en el que se contemplan procesos
# de difusion y reacciones quimicas.


def actualizacion_as(a, X, d):
    for i in range(1,d + 1):
        for j in range(1,d + 1):
            a[0,i,j] =  D*(X[0,i,j] - X[0,i - 1,j])
            if a[0,i,j] < 0:
                a[0][i][j] = 0
            a[1,i,j] =  D*(X[0,i,j] - X[0,i + 1,j])
            if a[1,i,j] < 0:
                a[1,i,j] = 0
            a[2,i,j] =  D*(X[0,i,j] - X[0,i,j - 1])
            if a[2,i,j] < 0:
                a[2,i,j] = 0
            a[3,i,j] =  D*(X[0,i,j] - X[0,i,j + 1])
            if a[3,i,j] < 0:
                a[3,i,j] = 0
                
            # R_mu: S1 -> S2
            a[4,i,j] = K1*X[0,i,j]
            #R_nu: S2 -> S1
            a[5,i,j] = K2*X[1,i,j]
    return a

# Funcion que suma todas las propensiones para poder calcular las probabilidades
# Probabilidades por unidad de tiempo


def fun_a0(a, d, cr):
    a0 = 0
    for mu in range(0, cr):
        for i in range(1,d + 1):
            for j in range(1,d + 1):
                a0 = a0 + a[mu][i][j]
    return a0

# Funcion que actualiza las condiciones de frontera dado un vector X concentraciones para
# todas las especies quimicas presentes en X.


def func_act_CC(X0, d):
    n_chs = X0.shape[0]
    for s in range(n_chs):
        for j in range(1,d + 1):
            X0[s][0][j] = X0[s][d][j]
            X0[s][d + 1][j] = X0[s][1][j]
            X0[s][j][0] = X0[s][j][d]
            X0[s][j][d + 1] = X0[s][j][1]
    return X0

# Funcion que actualiza las probabilidades usando el vector de propensiones a
# como parametros d, el tamaño de la malla, cr, el numero de reacciones o procesos
# de difusion

def act_probs(a, a0, P2, P3, d, cr):

    for mu in range(cr):
        for i in range(1,d + 1):
            for j in range(1,d + 1):
                P2[mu][i][j] = a[mu][i][j]/a0
                
    for mu in range(cr):
        suma_mu = 0
        for i in range(1,d + 1):
            for j in range(1,d + 1):
                suma_mu = suma_mu + P2[mu][i][j]
        P3[mu] = suma_mu

    return P2, P3

# Funcion que elige un proceso de reaccion o difusion
# utilizando P3 y un numero pseudoaleatorio

def fun_mu(P3,r2,cr):
    suma_mu = 0.
    for nu in range(cr):
        #print('P3[mu]', P3[nu],'nu',nu)
        suma_mu = suma_mu + P3[nu]
        #print('Suma_mu', suma_mu,'r2',r2)
        if suma_mu >= r2:
            return nu
    return nu

# Funcion que elige una vez que se haya elegido un proceso de reaccion
# o de difusion un renglon donde ocurre

def fun_i(P3,P4,r3,mu,d):
    suma_i = 0
    for i_star in range(1,d + 1):
        suma_i = suma_i + P4[i_star]
        if suma_i >= r3*P3[mu]:
            return i_star
    return i_star

# Funcion que elige, una vez elegido el proceso de reaccion o difusion
# y un renglon, una columna de la malla


def fun_j(P4,r4,mu,i,d):
    suma_j = 0.
    for j_star in range(1,d + 1):
        suma_j = suma_j + P2[mu][i][j_star]
        if suma_j >= r4*(P4[i]):
            return j_star
    return j_star


# Esta es otra de las funciones medulares para el algoritmo de Gillespie
# Actualiza las poblaciones moleculares de acuerdo al proceso ocurrido y
# al lugar de ocurrencia

def fun_act_dif(X, mu, i, j, d):
    # Difusion hacia arriba
    if mu == 0:
        X[0,i,j] = X[0,i,j] - 1
        X[0,i - 1, j] = X[0,i - 1, j] + 1
        if i == 1:
            X[0,d, j] = X[0,d, j] + 1
    # Difusion hacia abajo    
    elif mu == 1:
        X[0,i,j] = X[0,i,j] - 1
        X[0,i + 1, j] = X[0,i + 1, j] + 1
        if i == d:
            X[0,1, j] = X[0,1, j] + 1 
    # Difusion hacia la izquierda
    elif mu == 2:
        X[0,i,j] = X[0,i,j] - 1
        X[0,i,j - 1] = X[0,i,j - 1] + 1
        if j == 1:
            X[0,i, d] = X[0,i, d] + 1
    # Difusion hacia la derecha
    elif mu == 3:
        X[0,i,j] = X[0,i,j] - 1
        X[0,i,j + 1] = X[0,i,j + 1] + 1
        if j == d:
            X[0,i,1] = X[0,i,1] + 1
    elif mu == 4:
        # R_mu: S1 -> S2
        X[0,i,j] = X[0,i,j] - 1
        X[1,i,j] = X[1,i,j] + 1
    elif mu == 5:
        # R_mu: S2 -> S1
        X[1,i,j] = X[1,i,j] - 1
        X[0,i,j] = X[0,i,j] + 1
        
    return X

# Funcion que integra todas las sub rutinas de un paso montecarlo
# Avanza el tiempo, actualiza las poblaciones moleculares, las propensiones
# y las probabilidades en cada paso


def paso_montecarlo(X, a, t, P2, P3, d, cr):
    
    r1 = np.random.rand()
    r2 = np.random.rand()
    r3 = np.random.rand()
    r4 = np.random.rand()
    
    
    # Calculo del peso universal a0
    
    a0 = fun_a0(a, d, cr)
    
    #Calculo del intervalo tau
    
    
    tau = (1./a0)*np.log(1/r1)
    
    #Actualizamos el tiempo
                
    t = t + tau
    
    
    # Seleccion del procesos de reaccion-difusion
    
    mu = fun_mu(P3,r2,cr)
    
        # Aqui se tienen que calcular las probabilidades P4 y P5
    
    for i in range(1,d + 1):
        suma_j = 0
        for j in range(1,d + 1):
            suma_j = suma_j + P2[mu,i,j]
        P4[i] = suma_j
        
    # Seleccion del renglon de ocurrencia
        
    i = fun_i(P3,P4,r3,mu,d)
    
    
    
    # Seleccion de la columna de ocurrencia
    
    j = fun_j(P4,r4,mu,i,d)
    
    #print('tau',format(tau, '.4f'),'mu',mu,'i',i,'j',j)
    
    #print('X[i,j]',X[i][j], 'X[i - 1,j]',X[i - 1][j], 'X[i + 1,j]',X[i + 1][j],#
    #      'X[i,j - 1]',X[i][j - 1], 'X[i,j + 1]',X[i][j + 1], 'a(mu,i,j)',a[mu][i][j],#
    #      'P2[mu][i][j]',P2[mu][i][j],'P3[mu]',P3[mu], 'P4[mu][i]',P4[i])
    
    # Actualizamos el vector de concentracion de acuerdo con el proceso elegido
    
    #print(a[4,5,5],a[5,5,5])

    X = fun_act_dif(X, mu, i, j, d)
    
    #print('X[i,j]',X[i][j])
    
    # Actualizamos las condiciones de contorno si se actualizo algun punto de la frontera
    
    #if i == 1 | i == d | j == 1 | j == d:
    X = func_act_CC(X,d)

        
    # Actualizamos las as de acuerdo al nuevo vector X de concentracion
    
    a = actualizacion_as(a, X, d)
    
    # Actualizamos el peso universal a0
    
    a0 = fun_a0(a, d, cr)
    
    
    # Actualizamos las Ps
    
    P2, P3 = act_probs(a,a0,P2,P3, d,cr)
    

    
    return t, X, a, P2, P3

