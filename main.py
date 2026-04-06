from Shamir import Shamir, ShamirSimplificado
from Rampa import ShamirRampa, McElieceSarwate
from galois import GF

from typing import Callable

def pedir_entero(pregunta: str, mensaje_error: str, condicion: Callable[[int], bool]) -> int:
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

cuerpo = GF(2 ** 64)

# Pedir los participantes y los datos del esquema
participantes = []
n = pedir_entero('Escriba el número de particiantes: ',
                   f'El número de participantes debe ser al menos 2.', lambda x: 2 <= x)
for i in range(n):
    participante = input(f'Escriba el nombre del participante nº{i + 1}: ')
    while participante in participantes:
        participante = input(f'El participante introducido ya existe, por favor introduzca otro: ')
    participantes.append(participante)
print() # Imprimir espacio en blaco para mejor claridad visual

r = pedir_entero('Escriba número de participantes necesarios para recuperar el secreto: ',
                   f'número de participantes necesarios para recuperar el secreto debe ser al menos 2 y menor o igual que el número de participantes ({n})', lambda x: 2 <= x <= n )

l = pedir_entero('Escriba el número de secretos que se desea compartir: ',
                   f'El número de secretos debe ser al menos 1 y menor que el parámetro de reconstrucción ({r}).', lambda x: 1 <= x < r)

# Preguntar por el esquema que se desea usar
if l == 1:
    if r == n:
        yn = input('Se ha detectado que el número de participantes necesarios para reconstruir el esquema coincide con el número de participantes, ¿desea utilizar el esquema simplificado? (y/n): ')
        if yn.lower() in ('si', 's', 'y', 'yes'):
            ss = ShamirSimplificado(cuerpo, participantes)
        else:
            ss = Shamir(cuerpo, r, participantes)
    else:
        ss = Shamir(cuerpo, r, participantes)
else:
    print('¿Cual de los dos siguientes esquemas desea usar?\n1.- Esquema de Shamir en rampa.\n2.- Esquema de McEliece-Sarwate.')
    esq = pedir_entero('Respuesta: ',
                 f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
    if esq == 1:
        ss = ShamirRampa(cuerpo, r, l, participantes)
    else:
        ss = McElieceSarwate(cuerpo, r, l, participantes)

# Preguntar por participaciones anticipadas
yn = input('¿Desea repartir participación anticipada? (y/n): ')
if yn.lower() in ('si', 's', 'y', 'yes'):
    n_anticipados = pedir_entero(f'Introduzca el número de participaciones anticipadas (1 - {r-l}): ',
                 f'El número de participaciones anticipadas debe ser al menos 1 y menor que el parámetro de privacidad {r-l}', lambda x: 1 <= x <= r-l)
    participantes_anticipados = []
    # Añadir participantes anticipados
    for i in range(n_anticipados):
        part_anticipado = input(f'Escriba el nombre del participante anticipado nº{i + 1}: ')
        # Comporbar validez de los participantes anticipados
        while True:
            if part_anticipado in participantes_anticipados:
                print('Se ha introducido un participante duplicado.')
            elif part_anticipado not in participantes:
                print('Se ha introducido un participante no registrado.')
            else:
                break
            part_anticipado = input(f'Escriba el nombre del participante anticipado nº{i + 1}: ')
        participantes_anticipados.append(part_anticipado)
    # Crear particpaciones anticipadas
    participaciones_anticipadas = ss.crear_anticipadas(participantes_anticipados)
    for nombre, ant_b64 in participaciones_anticipadas:
        print(f'{nombre}: {ant_b64}')

# Pedir el secreto y codificarlo
if l == 1:
    secreto = input(f'Escriba el secreto: ').encode()
    while len(secreto) > 8:
        secreto = input('Se ha introducido un secreto con longitud de bytes mayor que 8, introduzca otro: ').encode()
else:
    secreto = []
    for i in range(l):
        sec = input(f'Escriba el secreto nº{i + 1}: ').encode()
        while len(sec) > 8:
            sec = input(
                'Se ha introducido un secreto con longitud de bytes mayor que 8, introduzca otro: ').encode()
        secreto.append(sec)

# Crear participaciones
participaciones = ss.crear_participaciones(secreto)
for nombre, part_b64 in participaciones:
    print(f'{nombre}: {part_b64}')

# Reconstrucción del secreto
print('Escriba el nombre de los participantes que busca reconstruir el secreto y su participacion.')
conjunto = []
nombres = []
# Pedir las participaciones
for _ in range(r):
    nombre = input('Nombre: ')
    # Comprobar validez de participantes
    while True:
        if nombre in nombres:
            print('Se ha introducido un participante duplicado.')
        elif nombre not in participantes:
            print('Se ha introducido un participante no registrado.')
        else:
            break
        nombre = input('Nombre: ')
    participacion = input('Participacion: ')
    conjunto.append((nombre, participacion))
    nombres.append(nombre)

# Reconstruir el secreto
secreto = ss.recuperar_secreto(conjunto)

try:
    if l == 1:
        print('secreto:', secreto.decode())
    else:
        print('secreto:', [s.decode() for s in secreto])
except UnicodeDecodeError:
    print('El secreto obtenido no es válido')