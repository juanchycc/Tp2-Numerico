import numpy as np

def leer_archivo():
	cantidad_elementos = 0
	temperatura = []
	variacion_iodo = []

	with open("data_00.txt","r") as file:                         
		for linea in file.readlines():      
			linea = linea.split()
			temperatura.append(float(linea[0]))
			variacion_iodo.append(float(linea[1]))
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

#por ahora no la use era mejor hacer todos los prod vactoriales en la misma iteracion
def producto_vectorial(primer_vector,segundo_vector,cantidad_elementos):
	aux = 0
	for i in range(cantidad_elementos):
		aux += (primer_vector[i]*segundo_vector[i])
	return aux



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

	return error

#linealiza como ln( 1/SI - 1 ) = ( ln(a) . 1 ) - (b . T) --> no funciona, por el logaritmo con negativo
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


#item c
def ajuste_polinomial_tres(temperatura,variacion_iodo,cantidad_elementos):
	fi_cero = []
	fi_uno = temperatura
	fi_dos = []
	fi_tres = []
	H = variacion_iodo
	filas = 4
	columnas = 5

	for i in range(cantidad_elementos):
		fi_cero.append(1)
		fi_dos.append( ( fi_uno[i] * fi_uno[i] ) )
		fi_tres.append( ( fi_dos[i] * fi_uno[i] ) )

	matriz = inicializar_matriz(filas,columnas)

	for i in range(cantidad_elementos):
		matriz[0][0] += ( fi_cero[i] * fi_cero[i] )
		matriz[0][1] += ( fi_cero[i] * fi_uno[i] )
		matriz[0][2] += ( fi_cero[i] * fi_dos[i] )
		matriz[0][3] += ( fi_cero[i] * fi_tres[i] )
		matriz[0][4] += ( fi_cero[i] * H[i] )

		matriz[1][0] =	matriz[0][1]
		matriz[1][1] += ( fi_uno[i] * fi_uno[i] )
		matriz[1][2] += ( fi_uno[i] * fi_dos[i] )
		matriz[1][3] += ( fi_uno[i] * fi_tres[i] ) 
		matriz[1][4] += ( fi_uno[i] * H[i] )

		matriz[2][0] = matriz[0][2]
		matriz[2][1] = matriz[1][2]
		matriz[2][2] += ( fi_dos[i] * fi_dos[i] )
		matriz[2][3] += ( fi_dos[i] * fi_tres[i] )
		matriz[2][4] += ( fi_dos[i] * H[i] )

		matriz[3][0] = matriz[0][3]
		matriz[3][1] = matriz[1][3]
		matriz[3][2] = matriz[2][3]
		matriz[3][3] += ( fi_tres[i] * fi_tres[i] )
		matriz[3][4] += ( fi_tres[i] * H[i] )

	#método de Crout (descomposición LU)
	#triangulación
	m21 = ( matriz[1][0] / matriz[0][0] )
	m31 = ( matriz[2][0] / matriz[0][0] )
	m41 = ( matriz[3][0] / matriz[0][0] )

	for i in range(columnas):
		matriz[1][i] = ( ( matriz[0][i] * m21 ) - matriz[1][i] )
		matriz[2][i] = ( ( matriz[0][i] * m31 ) - matriz[2][i] )
		matriz[3][i] = ( ( matriz[0][i] * m41 ) - matriz[3][i] )
		
	m32 = ( matriz[2][1] / matriz[1][1] )
	m42 = ( matriz[3][1] / matriz[1][1] )

	x = 1
	for x in range(columnas):
		matriz[2][x] = ( ( matriz[1][x] * m32 ) - matriz[2][x] )
		matriz[3][x] = ( ( matriz[1][x] * m42 ) - matriz[3][x] )

	m43 = ( matriz[3][2] / matriz[2][2])
	x = 2
	for x in range(columnas):
		matriz[3][x] = ( ( matriz[2][x] * m43 ) - matriz[3][x] )

	L = inicializar_matriz(4,4)
	U = matriz

	L[0][0] = 1
	L[1][0] = m21
	L[1][1] = 1
	L[2][0] = m31
	L[2][1] = m32
	L[2][2] = 1
	L[3][0] = m41
	L[3][1] = m42
	L[3][2] = m43
	L[4][2] = 1

	#  L . Y = b --> L . Y = H
	Y = []
	Y[0] = ( H[0] / L[0][0] )
	Y[1] = ( ( H[1] - ( L[1][0] * Y[0] ) ) / L[1][1] )
	Y[2] = ( ( H[2] - ( L[2][0] * Y[0]) - ( L[2][1] * Y[1] ) ) / L[2][2] )
	Y[3] = ( ( H[3] - ( L[3][0] * Y[0]) - ( L[3][1] * Y[1]) - (L[3][2] * Y[2]) ) / L[3][3] )

	# U . X = Y
	X = []
	X[3] = ( Y[3] / U[3][3] )
	X[2] = ( ( Y[2] - ( U[2][3] * X[3] ) ) / U[2][2] )
	X[1] = ( ( Y[1] - ( U[1][3] * X[3] ) - ( U[1][2] * X[2] ) ) / U[1][1] )
	X[0] = ( ( Y[0] - ( U[0][3] * X[3] ) -( U[0][2] * X[2] ) - ( U[0][1] * X[1] ) ) / U[0][0])
	

def main():
	temperatura,variacion_iodo,cantidad_elementos = leer_archivo()

	ajuste_lineal(temperatura,variacion_iodo,cantidad_elementos)
	ajuste_polinomial(temperatura,variacion_iodo,cantidad_elementos)
	#ajuste_logistico(temperatura,variacion_iodo,cantidad_elementos)
	ajuste_polinomial_tres(temperatura,variacion_iodo,cantidad_elementos)
main()