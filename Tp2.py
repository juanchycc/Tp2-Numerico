

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

def main():
	temperatura,variacion_iodo,cantidad_elementos = leer_archivo()

	ajuste_lineal(temperatura,variacion_iodo,cantidad_elementos)
	ajuste_polinomial(temperatura,variacion_iodo,cantidad_elementos)
main()