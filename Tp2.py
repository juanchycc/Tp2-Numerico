from os import X_OK
import numpy as np
import matplotlib.pyplot as plt

#funciones auxiliares
def leer_archivo(nombre_archivo):
	cantidad_elementos = 0
	temperatura = []
	variacion_iodo = []

	with open(nombre_archivo,"r") as file:                         
		for linea in file.readlines():      
			linea = linea.split()
			temperatura.append(float(linea[0]))
			variacion_iodo.append(float(linea[1])) #cambio de unidad (lo saque porque fallecia incluso antes)
			cantidad_elementos += 1
		

	return temperatura,variacion_iodo,cantidad_elementos

#recibe el tamaño de fila y columna y devuelve una matriz de esa dimension inicializada en ceros
def inicializar_matriz(filas,columnas):
	matriz = []
	for i in range(filas):
		matriz.append([])
		for j in range(columnas):
			matriz[i].append(0)
	return matriz

def printear_matriz(matriz): # es solo para printear ordenadamente
	for i in range(len(matriz)):
		for j in range(len(matriz[0])):
			print(matriz[i][j],end=" ")
		print()

#item a
def ajuste_lineal(temperatura,variacion_iodo,cantidad_elementos):
	fi_cero = []
	fi_uno = temperatura
	H = variacion_iodo
	filas = 2
	columnas = 3

	for i in range(cantidad_elementos):
		fi_cero.append(1)

	matriz = inicializar_matriz(filas,columnas)

	for i in range(cantidad_elementos):
		matriz[0][0] += ( fi_cero[i] * fi_cero[i] )
		matriz[0][1] += ( fi_cero[i] * fi_uno[i] )
		matriz[0][2] += ( fi_cero[i] * H[i] )
		matriz[1][0] =	matriz[0][1]
		matriz[1][1] += ( fi_uno[i] * fi_uno[i] )
		matriz[1][2] += ( fi_uno[i] * H[i] )
	
	#eliminacion de gauss en fc externa?
	m21 = ( matriz[1][0] / matriz[0][0] )
	for i in range(columnas):
		matriz[1][i] = ( ( matriz[0][i] * m21 ) - matriz[1][i] )

	b = ( matriz[1][2] / matriz[1][1] )
	a = ( ( matriz[0][2] - ( b * matriz[0][1])) / matriz[0][0] )
	#calculo el error
	error = 0
	for i in range(cantidad_elementos):
		aux = ( variacion_iodo[i] - ( a + ( b * temperatura[i] ) ) )
		error += ( aux * aux )

	print(error)
	return error

def ajuste_polinomial(temperatura,variacion_iodo,cantidad_elementos):
	fi_cero = []
	fi_uno = temperatura
	fi_dos = []
	H = variacion_iodo
	filas = 3
	columnas = 4

	for i in range(cantidad_elementos):
		fi_cero.append(1)
		fi_dos.append( (fi_uno[i] * fi_uno[i]))

	matriz = inicializar_matriz(filas,columnas)

	for i in range(cantidad_elementos):
		matriz[0][0] += ( fi_cero[i] * fi_cero[i] )
		matriz[0][1] += ( fi_cero[i] * fi_uno[i] )
		matriz[0][2] += ( fi_cero[i] * fi_dos[i] )
		matriz[0][3] += ( fi_cero[i] * H[i] )

		matriz[1][0] =	matriz[0][1]
		matriz[1][1] += ( fi_uno[i] * fi_uno[i] )
		matriz[1][2] += ( fi_uno[i] * fi_dos[i] ) 
		matriz[1][3] += ( fi_uno[i] * H[i] )

		matriz[2][0] = matriz[0][2]
		matriz[2][1] = matriz[1][2]
		matriz[2][2] += ( fi_dos[i] * fi_dos[i] )
		matriz[2][3] += ( fi_dos[i] * H[i] )

	#gauss
	m21 = ( matriz[1][0] / matriz[0][0] )
	m31 = ( matriz[2][0] / matriz[0][0] )

	for i in range(columnas):
		matriz[1][i] = ( ( matriz[0][i] * m21 ) - matriz[1][i] )
		matriz[2][i] = ( ( matriz[0][i] * m31 ) - matriz[2][i] )

	m32 = ( matriz[2][1] / matriz[1][1] )

	x = 1
	for x in range(columnas):
		matriz[2][x] = ( ( matriz[1][x] * m32 ) - matriz[2][x] )

	c = ( matriz[2][3] / matriz[2][2] )
	b = ( ( matriz[1][3] - ( matriz[1][2] * c ) ) / matriz[1][1] )
	a = ( ( matriz[0][3] - ( matriz[0][2] * c ) - ( matriz[0][1] * b ) ) / matriz[0][0] )

	#error
	error = 0
	for i in range(cantidad_elementos):
		aux = ( variacion_iodo[i] - ( a + ( b * temperatura[i] ) + ( c * fi_dos[i] ) ) )
		error += ( aux * aux )

	print(error)	
	return error

#linealiza como ln( 1/SI - 1 ) = ( ln(a) . 1 ) - (b . T) 
def ajuste_logistico(temperatura,variacion_iodo,cantidad_elementos):
	fi_cero = []
	fi_uno = temperatura
	H = []
	filas = 2
	columnas = 3
	matriz = inicializar_matriz(filas,columnas)

	for i in range(cantidad_elementos):
		H.append(np.log(  (1/variacion_iodo[i]) - 1 ))
		fi_cero.append(1)

	#gauss idem a ajuste lineal, por matriz de 2x3
	for i in range(cantidad_elementos):
		matriz[0][0] += ( fi_cero[i] * fi_cero[i] )
		matriz[0][1] += ( fi_cero[i] * fi_uno[i] )
		matriz[0][2] += ( fi_cero[i] * H[i] )
		matriz[1][0] =	matriz[0][1]
		matriz[1][1] += ( fi_uno[i] * fi_uno[i] )
		matriz[1][2] += ( fi_uno[i] * H[i] )

	m21 = ( matriz[1][0] / matriz[0][0] )
	for i in range(columnas):
		matriz[1][i] = ( ( matriz[0][i] * m21 ) - matriz[1][i] )

	b = ( matriz[1][2] / matriz[1][1] )
	c0 = ( ( matriz[0][2] - ( b * matriz[0][1])) / matriz[0][0] )

	# ln(a) = c0 --> a = exp(c0)
	a =  np.e**c0
	error = 0
	for i in range(cantidad_elementos):
		aux = ( variacion_iodo[i] - ( 1 / ( 1 + ( a * ( np.e**(-b*temperatura[i])))) ) )
		error += ( aux * aux )

	print(error)
	return error

#item b
def ajuste_cuadratico_sin_linealizar(temperatura,variacion_iodo,cantidad_elementos,tolerancia):
	solucion2 = np.matrix([[190],[10],[0.05]])
	error = 10000
	contador = 0
	while error > tolerancia:
		a2 = float(solucion2[0][0])		
		b2 = float(solucion2[1][0])
		c2 = float(solucion2[2][0])
		
		valores_funcion2 = funcion_sigmoidea(a2,b2,c2,temperatura,variacion_iodo,cantidad_elementos)
		jacobiano2 = obtener_jacobiano2(a2,b2,c2,temperatura,variacion_iodo,cantidad_elementos)
		jacobiano2 = np.matrix(jacobiano2)
		valores_funcion2 = np.matrix([[valores_funcion2[0]],[valores_funcion2[1]],[valores_funcion2[2]]])
		delta2 = (jacobiano2**-1)*(valores_funcion2*-1)
		vieja_solucion_2 = solucion2		
		solucion2 = vieja_solucion_2 + delta2
		
		a = float(solucion2[0][0]) 
		b = float(solucion2[1][0])
		c = float(solucion2[2][0])

		error = np.sqrt((a-a2)**2 + (b-b2)**2 + (c-c2)**2)	
		contador += 1

		print(contador)
	return solucion2

def funcion_sigmoidea(a,b,c,temperatura,variacion_iodo,cantidad_elementos):
	valores_funcion = [0,0,0]
	for i in range(cantidad_elementos):
		x = temperatura[i]
		y = variacion_iodo[i]
		valores_funcion[0] += 2*(y- a*np.e**(-b*np.e**(-c*x))) *  np.e**(-b*np.e**(-c*x)) *(-1)
		valores_funcion[1] += 2*(y - a*np.e**(-b*np.e**(-c*x))) *  a*np.e**(-c*x-b*np.e**(-c*x)) #*(-1) aca corregi cosas tambien
		valores_funcion[2] += 2*(y - a*np.e**(-b*np.e**(-c*x))) *	-a*b*x*np.e**(-b*np.e**(-c*x)-c*x)  #error en esta fórmula
	
	return valores_funcion

def obtener_jacobiano2(a,b,c,temperatura,variacion_iodo,cantidad_elementos):
	jacobiano = [[0,0,0],[0,0,0],[0,0,0]]
	for i in range(cantidad_elementos):
		x = temperatura[i]
		d = variacion_iodo[i]                

		jacobiano [0][0] +=  2*np.e**(-2*b*np.e**(-c*x))  #1
		jacobiano [0][1] += -2*(-np.e**(-np.e**(-c*x)*b-c*x)*d+2*a*np.e**(-2*np.e**(-c*x)*b-c*x))  #2
		jacobiano [0][2] += -2*(-2*a*b*np.e**(-2*b*np.e**(-c*x)-c*x)*x+b*np.e**(-b*np.e**(-c*x)-c*x)*d*x) #3
		jacobiano [1][0] +=  2*np.e**(-b*np.e**(-c*x)-c*x)*(d-2*np.e**(-b*np.e**(-c*x))*a)   #4
		jacobiano [1][1] +=  2*a*(-np.e**(-np.e**(-c*x)*b-2*c*x)*d+2*a*np.e**(-2*np.e**(-c*x)*b-2*c*x))
		jacobiano [1][2] +=  2*a*(a*np.e**(-2*b*np.e**(-c*x)-c*x)*x-2*a*b*np.e**(-2*b*np.e**(-c*x)-2*c*x)*x+b*np.e**(-b*np.e**(-c*x)-2*c*x)*d*x-np.e**(-b*np.e**(-c*x)-c*x)*d*x)
		jacobiano [2][0] += -2*b*np.e**(-b*np.e**(-c*x)-c*x)*x*(-2*np.e**(-b*np.e**(-c*x))*a+d)
		jacobiano [2][1] += -2*a*x*(2*a*np.e**(-2*np.e**(-c*x)*b-2*c*x)*b-np.e**(-np.e**(-c*x)*b-2*c*x)*d*b-a*np.e**(-2*np.e**(-c*x)*b-c*x)+np.e**(-np.e**(-c*x)*b-c*x)*d)		#8 nuevo
		jacobiano [2][2] += -2*a*b*x*(a*np.e**(-2*b*np.e**(-c*x)-c*x)*x-2*a*b*np.e**(-2*b*np.e**(-c*x)-2*c*x)*x+b*np.e**(-b*np.e**(-c*x)-2*c*x)*d*x-np.e**(-b*np.e**(-c*x)-c*x)*d*x)  #9 nuevo

	return jacobiano

#item c
def ajuste_polinomial_tres(temperatura,variacion_iodo,cantidad_elementos):
	fi_cero = []
	fi_uno = temperatura
	fi_dos = []
	fi_tres = []
	H = variacion_iodo
	filas = 4
	columnas = 4

	for i in range(cantidad_elementos):
		fi_cero.append(1)
		fi_dos.append( ( fi_uno[i] * fi_uno[i] ) )
		fi_tres.append( ( fi_dos[i] * fi_uno[i] ) )

	A = inicializar_matriz(filas,columnas)

	for i in range(cantidad_elementos):
		A[0][0] += ( fi_cero[i] * fi_cero[i] )
		A[0][1] += ( fi_cero[i] * fi_uno[i] )
		A[0][2] += ( fi_cero[i] * fi_dos[i] )
		A[0][3] += ( fi_cero[i] * fi_tres[i] )

		A[1][0] =	A[0][1]
		A[1][1] += ( fi_uno[i] * fi_uno[i] )
		A[1][2] += ( fi_uno[i] * fi_dos[i] )
		A[1][3] += ( fi_uno[i] * fi_tres[i] ) 

		A[2][0] = A[0][2]
		A[2][1] = A[1][2]
		A[2][2] += ( fi_dos[i] * fi_dos[i] )
		A[2][3] += ( fi_dos[i] * fi_tres[i] )

		A[3][0] = A[0][3]
		A[3][1] = A[1][3]
		A[3][2] = A[2][3]
		A[3][3] += ( fi_tres[i] * fi_tres[i] )

	aux1 = inicializar_matriz(filas,columnas)
	aux2 = inicializar_matriz(filas,columnas)
	for i in range(filas):
		for j in range(columnas):
			aux1[i][j] = A[i][j]
			aux2[i][j] = A[i][j]

	print(factorizacion_doolitle(aux1,H,3,columnas))
	print(factorizacion_crout(A,H,3,columnas))
	print(factorizacion_doolitle(aux2,H,6,columnas))
	print(factorizacion_crout(A,H,6,columnas))

def factorizacion_crout(A,B,cantidad_digitos,columna):
	#resuelvo las ecuaciones para obtener L y U	
	for i in range(columna):
		B[i] = round(B[i],cantidad_digitos)
		for j in range(columna):
			A[i][j] = round(A[i][j],cantidad_digitos)

	L = inicializar_matriz(4,4)
	U = inicializar_matriz(4,4)

	U[0][0] = 1
	U[1][1] = 1
	U[2][2] = 1
	U[3][3] = 1
	L[0][0] = A[0][0]
	U[0][1] = round(( A[0][1] / L[0][0] ),cantidad_digitos)
	U[0][2] = round(( A[0][2] / L[0][0] ),cantidad_digitos)
	U[0][3] = round(( A[0][3] / L[0][0] ),cantidad_digitos)
	L[1][0] = A[1][0]
	L[1][1] = round((A[1][1] - round((L[1][0] * U[0][1]),cantidad_digitos)),cantidad_digitos)
	U[1][2] = round((round((A[1][2] - round((L[1][0] * U[0][1]),cantidad_digitos)),cantidad_digitos) / L[1][1]),cantidad_digitos)
	U[1][3] = round(((round((A[2][1] - round((L[2][0] * U[0][1]),cantidad_digitos)),cantidad_digitos)) / L[1][0]),cantidad_digitos)
	L[2][0] = A[2][0]
	L[2][1] = round((A[2][1] - round((L[2][0] * U[0][1]),cantidad_digitos)),cantidad_digitos)
	L[2][2] = round((A[2][2] - round((L[2][0] * U[0][1]),cantidad_digitos) - round((L[2][1] * U[1][2]),cantidad_digitos)),cantidad_digitos)
	U[2][3] = round(((round((A[2][2] - round((L[2][0] * U[0][3]),cantidad_digitos) - round((L[2][1] * U[1][3]),cantidad_digitos)),cantidad_digitos)) / L[2][2]),cantidad_digitos)
	L[3][0] = A[3][0]
	L[3][1] = round((A[3][1] - round((L[3][0] * U[0][1]),cantidad_digitos)),cantidad_digitos)
	L[3][2] = round((A[3][2] - round((L[3][0] * U[0][2]),cantidad_digitos) - round((L[3][1] * U[1][2]),cantidad_digitos)),cantidad_digitos)
	L[3][3] = round((A[3][3] - round((L[3][0] * U[0][3]),cantidad_digitos) - round((L[3][1] * U[1][3]),cantidad_digitos) - round((L[3][2] * U[2][3]),cantidad_digitos)),cantidad_digitos)

	return resolver_sistema(L,U,B,cantidad_digitos)

def factorizacion_doolitle(A,B,cantidad_digitos,columna):
	for i in range(columna):
		B[i] = round(B[i],cantidad_digitos)
		for j in range(columna):
			A[i][j] = round(A[i][j],cantidad_digitos)

	m21 = round(( A[1][0] / A[0][0] ),cantidad_digitos)
	m31 = round(( A[2][0] / A[0][0] ),cantidad_digitos)
	m41 = round(( A[3][0] / A[0][0] ),cantidad_digitos)

	for i in range(columna):
		A[1][i] = round(( round(( A[0][i] * m21 ),cantidad_digitos) - A[1][i] ),cantidad_digitos)
		A[2][i] = round(( round(( A[0][i] * m31 ),cantidad_digitos) - A[2][i] ),cantidad_digitos)
		A[3][i] = round(( round(( A[0][i] * m41 ),cantidad_digitos) - A[3][i] ),cantidad_digitos)
		
	m32 = round(( A[2][1] / A[1][1] ),cantidad_digitos)
	m42 = round(( A[3][1] / A[1][1] ),cantidad_digitos)

	x = 1
	for x in range(columna):
		A[2][x] = round(( round(( A[1][x] * m32 ),cantidad_digitos) - A[2][x] ),cantidad_digitos)
		A[3][x] = round(( round(( A[1][x] * m42 ),cantidad_digitos) - A[3][x] ),cantidad_digitos)

	m43 = round(( A[3][2] / A[2][2]),cantidad_digitos)
	x = 2
	for x in range(columna):
		A[3][x] = round(( round(( A[2][x] * m43 ),cantidad_digitos) - A[3][x] ),cantidad_digitos)

	L = inicializar_matriz(4,4)
	U = A

	L[0][0] = 1
	L[1][0] = m21
	L[1][1] = 1
	L[2][0] = m31
	L[2][1] = m32
	L[2][2] = 1
	L[3][0] = m41
	L[3][1] = m42
	L[3][2] = m43
	L[3][3] = 1
	return resolver_sistema(L,U,B,cantidad_digitos)

def resolver_sistema(L,U,B,cantidad_digitos):
	#  L . Y = b
	Y = [0,0,0,0]
	Y[0] = round(( B[0] / L[0][0] ),cantidad_digitos)
	Y[1] = round(( round(( B[1] - round(( L[1][0] * Y[0]),cantidad_digitos) ),cantidad_digitos) / L[1][1] ),cantidad_digitos)
	Y[2] = round(( round(( B[2] - round(( L[2][0] * Y[0]),cantidad_digitos) )- round(( L[2][1] * Y[1]),cantidad_digitos) ,cantidad_digitos) / L[2][2] ),cantidad_digitos)
	Y[3] = round(( round(( B[3] - round(( L[3][0] * Y[0]),cantidad_digitos) )- round(( L[3][1] * Y[1]),cantidad_digitos) - round((L[3][2] * Y[2]),cantidad_digitos ),cantidad_digitos) / L[3][3] ),cantidad_digitos)
	# U . X = Y
	X = [0,0,0,0]
	X[3] = round(( Y[3] / U[3][3] ),cantidad_digitos)
	X[2] = round(( round(( Y[2] - round(( U[2][3] * X[3] ),cantidad_digitos) ),cantidad_digitos) / U[2][2] ),cantidad_digitos)
	X[1] = round(( round(( Y[1] - round(( U[1][3] * X[3] ),cantidad_digitos) - round(( U[1][2] * X[2] ),cantidad_digitos) ),cantidad_digitos) / U[1][1] ),cantidad_digitos)
	X[0] = round(( round(( Y[0] - round(( U[0][3] * X[3] ),cantidad_digitos) - round(( U[0][2] * X[2] ),cantidad_digitos) - round(( U[0][1] * X[1] ),cantidad_digitos) ),cantidad_digitos) / U[0][0]),cantidad_digitos)
	return X


def graficar(sigmoidea,temperatura,variacion_iodo):
	a = float(sigmoidea[0][0]) 
	b = float(sigmoidea[1][0])
	c = float(sigmoidea[2][0])
	plt.figure(dpi = 125)
	plt.scatter(temperatura,variacion_iodo)
	temperatura = np.arange(20,60,0.1)
	plt.plot(temperatura,a*np.e**(-b*np.e**(-c*temperatura)))
	plt.show()


def main():

	tolerancia = 0.00000001 # (item b)
	temperatura,variacion_iodo,cantidad_elementos = leer_archivo("data_00.txt")
	variacion_iodo_achicada = []
	for i in range(cantidad_elementos):
		variacion_iodo_achicada.append( variacion_iodo[i] / 1000 )

	ajuste_lineal(temperatura,variacion_iodo_achicada,cantidad_elementos)
	ajuste_polinomial(temperatura,variacion_iodo_achicada,cantidad_elementos)
	ajuste_logistico(temperatura,variacion_iodo_achicada,cantidad_elementos)
	parametros = ajuste_cuadratico_sin_linealizar(temperatura,variacion_iodo,cantidad_elementos,tolerancia)
	ajuste_polinomial_tres(temperatura,variacion_iodo_achicada,cantidad_elementos)
	graficar(parametros,temperatura,variacion_iodo)
	



	temperatura,variacion_iodo,cantidad_elementos = leer_archivo("data_01.txt")
	variacion_iodo_achicada = []
	for i in range(cantidad_elementos):
		variacion_iodo_achicada.append( variacion_iodo[i] / 1000 )

	ajuste_lineal(temperatura,variacion_iodo_achicada,cantidad_elementos)
	ajuste_polinomial(temperatura,variacion_iodo_achicada,cantidad_elementos)
	ajuste_logistico(temperatura,variacion_iodo_achicada,cantidad_elementos)
	ajuste_cuadratico_sin_linealizar(temperatura,variacion_iodo,cantidad_elementos,tolerancia)
	ajuste_polinomial_tres(temperatura,variacion_iodo_achicada,cantidad_elementos)
	
main()
