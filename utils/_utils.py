import secrets
import numpy as np
from base64 import b64decode, b64encode
from galois import Poly

def bytes_to_int(string):
    """
    Obtiene la representacion en entero de unos bytes.
    Si se proporciona una secuencia de bytes, devuelve una lista de las representaciones de cada cadena.
    :param string: La cadena de caracteres cuya representacion en enetros se quiere obtener.
    :return: El numero entero que representa a la cadena proporcionada.
    """
    if isinstance(string, bytes):
        return int.from_bytes(string, byteorder='big')
    return list(int.from_bytes(b, byteorder='big') for b in string)

def int_to_bytes(numero):
    """
    Obtiene la representacion en bytes de un entero.
    Si se proporciona una secuancia de enteros, devuelve una lista de las representaciones de cada número.
    :param numero: El número entero cuya representación en string se quiere obtener.
    :return: La representacion en bytes del número entero.
    """
    if isinstance(numero, int):
        longitud = (numero.bit_length() + 7) // 8 if numero > 0 else 1
        return numero.to_bytes(longitud, byteorder='big')
    if isinstance(numero, np.ndarray):
        if numero.ndim == 0:
            n_int = int(numero)
            longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            return n_int.to_bytes(longitud, byteorder='big')
        representacion_str = []
        for n_int in numero.tolist():
            longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            representacion_str.append(n_int.to_bytes(longitud, byteorder='big'))
        return representacion_str
    representacion_str = []
    for n in numero:
        n_int = int(n)
        longitud = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
        representacion_str.append(n_int.to_bytes(longitud, byteorder='big'))
    return representacion_str

def int_to_b64str(lista_int, longitud):
    """
    Obtiene la codificación en base64 de cada número de un array de enteros.
    :param lista_int: El array de enteros.
    :param longitud: El número de bytes que se desean usar para representar cada número.
        Se crea un OverflowError si el número no se puede representar con el la longitud dada.
    :return: La lista con la codificación de cada entero.
    """
    if isinstance(lista_int, np.ndarray):
        return list(b64encode(numero.to_bytes(longitud, byteorder='big')).decode() for numero in lista_int.tolist())
    return list(b64encode(numero.to_bytes(longitud, byteorder='big')).decode() for numero in lista_int)

def b64str_to_int(lista_b64):
    """
    Obtiene la decodificación en base64 de cada string de una lista.
    :param lista_b64: La lista de strings codificados.
    :return: La lista con la decodificación de cada string.
    """
    return list(int.from_bytes(b64decode(str_b64), byteorder='big') for str_b64 in lista_b64)

def random_array(sup, n):
    """
    Devuelve un array de longitud n con números aleatorios criptograficamente seguros en el rango [0, sup).
    :param sup: Cota superior (no incluida) para los números aleatorios.
    :param n: La cabtidad de números deseado.
    :return: El array de números aleatorios.
    """
    return [secrets.randbelow(sup) for _ in range(n)]

def random_polinomial(cuerpo, grado):
    """
    Construye un polinomio aleatorio criptográficamente seguro sobre un cuerpo.
    :param cuerpo: Cuerpo sobre el que se construye el polinomio.
    :param grado: El grado del polinomio. Equivalentemente, El número de coeficientes del polinomio menos 1.
    :return: El polinomio aleatorio.
    """
    return Poly([secrets.randbelow(cuerpo.order) for _ in range(grado+1)], field=cuerpo)

def extend_matrix(matriz):
    """
    Obtiene la extensión de la matriz proporcionada formada por la representación matricial de cada uno de sus elementos en el cuerpo base y se convierte a formato booleano.
    :param matriz: Matriz sobre el cuerpo GF(2, m) a extender.
    :return: La matriz extendida sobre el cuerpo base, GF(2).
    """
    gf = type(matriz)  # Obtener el cuerpo de trabajo
    if gf.characteristic != 2:
        raise ValueError(f'Se esperaba un cuerpo con elemento base 2, pero se ha recibido {gf.characteristic}.')
    num_bits = gf.degree
    shape_final = np.array(matriz.shape) * num_bits # Obtener la forma final de la matriz
    base = gf.primitive_element ** np.arange(num_bits-1, -1, -1) # Base del cuerpo extendido, como vector es Big Endian, se toma la base Big Endian
    matriz_nueva_base = gf(np.kron(matriz, base)) # Pasar la matriz a la nueva base
    matriz_gf2 = matriz_nueva_base.vector() # Se toma la representacion de cada elemento de la matriz en el nuevo cuerpo
    matriz_gf2 = np.transpose(matriz_gf2, (0, 2, 1)).reshape(shape_final) # Alineación de dimensiones
    return matriz_gf2.view(np.ndarray).astype(bool) # Se devuelve la matriz extendida en formato booleano

def ask_int(pregunta, mensaje_error, condicion):
    """
    Pide por consola un número entero y valida su valor.
    :param pregunta: Mensaje que se muestra por pantalla para pedir el número.
    :param mensaje_error: Mensaje que se muestr apor pantalla si el número no cumple la condición.
    :param condicion: Función que recibe el entero introducido y devuelve True si es válido y False si no.
    :return: El número entero introducido y válido.
    """
    while True:
        try:
            num = int(input(pregunta))
        except ValueError:
            print('No se ha introducido un número válido.')
        else:  # Se ejecuta si no ha habido ninguna excepcion
            if not condicion(num):
                print(mensaje_error)
            else:
                print() # Imprimir espacio en blaco para mejor claridad visual
                return num