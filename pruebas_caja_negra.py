from esquemas_clasicos import Simplificado, Shamir, ShamirRampa, McElieceSarwate
from esquemas_cuanticos import CGL, Ogawa, ZhangMatsumoto

from qiskit.quantum_info import Statevector
from qiskit import QuantumRegister
from galois import GF

import numpy as np

nombres_esq = {
    'Simplificado': 'Esquema Simplificado',
    'Shamir': 'Esquema de Shamir',
    'ShamirRampa': 'Esquema de Shamir en rampa',
    'McElieceSarwate': 'Esquema de McEliece-Sarwate',
    'CGL': 'Esquema de Cleve-Gottesman-Lo',
    'Ogawa': 'Esquema de Ogawa et al.',
    'ZhangMatsumoto': 'Esquema de Zhang-Matsumoto',
}
clasico_umbral = [Shamir]
clasico_rampa = [ShamirRampa, McElieceSarwate]
cuantico_umbral = [CGL]
cuantico_rampa = [Ogawa, ZhangMatsumoto]


def validacion_datos(Esquema, cuerpo, r=None, l=None, participantes=None):
    print(f'[TEST] {nombres_esq[Esquema.__name__]}:')
    try:
        if Esquema in clasico_umbral + cuantico_umbral:
            Esquema(cuerpo, r, participantes)
        elif Esquema in clasico_rampa + cuantico_rampa:
            Esquema(cuerpo, r, l, participantes)
        else:  # Si no es ninguno de los anteriores es el esquema simplificado
            Simplificado(cuerpo, participantes)
    except Exception as e:
        print(f'{type(e).__name__}: {e}')
    print()

def validacion_secreto(Esquema, cuerpo, secreto, r=None, l=None, participantes=None):
    print(f'[TEST] {nombres_esq[Esquema.__name__]}:')
    if Esquema in clasico_umbral + cuantico_umbral:
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in clasico_rampa + cuantico_rampa:
        ss = Esquema(cuerpo, r, l, participantes)
    else:  # Si no es ninguno de los anteriores es el esquema simplificado
        ss = Simplificado(cuerpo, participantes)
    try:
        ss.codificacion(secreto)
    except Exception as e:
        print(f'{type(e).__name__}: {e}')
    print()

def validacion_participaciones(Esquema, cuerpo, r=None, l=None, participantes=None, opcion='ne'):
    print(f'[TEST] {nombres_esq[Esquema.__name__]}:')
    if Esquema in clasico_umbral:
        secreto = b'\x00'
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in cuantico_umbral:
        secreto = Statevector([1] + [0]*(cuerpo.order-1))
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in clasico_rampa:
        secreto = [b'\x00'] * l
        ss = Esquema(cuerpo, r, l, participantes)
    elif Esquema in cuantico_rampa:
        secreto = Statevector([1] + [0] * (cuerpo.order**l - 1))
        ss = Esquema(cuerpo, r, l, participantes)
    else:  # Si no es ninguno de los anteriores es el esquema simplificado
        secreto = b'\x00'
        ss = Simplificado(cuerpo, participantes)
    part = ss.codificacion(secreto)

    if opcion == 'ne': # Se entregan participaciones inexistentes
        if Esquema in cuantico_umbral + cuantico_rampa:
            extra = QuantumRegister(3, 'i')
        else:
            extra = ('i', 'AAAA')
        part = part[:r-1] + [extra]
    elif opcion == 'd': # Se entregan participaciones duplicados
        part = part[:r-1] + [part[0]]
    else: # Si no es ninguno, se entregan participaciones insuficientes
        part = part[:r-1]
    try:
        ss.decodificacion(part)
    except Exception as e:
        print(f'{type(e).__name__}: {e}')
    print()

def validacion_participaciones_anticipadas(Esquema, cuerpo, r=None, l=None, participantes=None, opcion='dp'):
    print(f'[TEST] {nombres_esq[Esquema.__name__]}:')
    if Esquema in clasico_umbral:
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in cuantico_umbral:
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in clasico_rampa:
        ss = Esquema(cuerpo, r, l, participantes)
    elif Esquema in cuantico_rampa:
        ss = Esquema(cuerpo, r, l, participantes)
    else:  # Si no es ninguno de los anteriores es el esquema simplificado
        ss = Simplificado(cuerpo, participantes)
    try:
        if opcion == 'ne': # Se entregan participantes inexistentes
            anticipados = participantes[:1] + ['i']
            ss.comparticion_anticipada(anticipados)
        elif opcion == 'd': # Se entregan participantes duplicados
            anticipados = participantes[:1] + participantes[:1]
            ss.comparticion_anticipada(anticipados)
        elif opcion == 'df': # Se entregan participantes duplicados en dos fases (solo clásicos)
            anticipados = participantes[:1]
            ss.comparticion_anticipada(anticipados)
            ss.comparticion_anticipada(anticipados)
        else: # Si no es ninguno, se entregan más de los esperados
            ss.comparticion_anticipada(participantes)
    except Exception as e:
        print(f'{type(e).__name__}: {e}')
    print()

def validacion_orden(Esquema, cuerpo, r=None, l=None, participantes=None, orden='dc'):
    print(f'[TEST] {nombres_esq[Esquema.__name__]}:')
    if Esquema in clasico_umbral:
        secreto = b'\x00'
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in cuantico_umbral:
        secreto = Statevector([1] + [0] * (cuerpo.order - 1))
        ss = Esquema(cuerpo, r, participantes)
    elif Esquema in clasico_rampa:
        secreto = [b'\x00'] * l
        ss = Esquema(cuerpo, r, l, participantes)
    elif Esquema in cuantico_rampa:
        secreto = Statevector([1] + [0] * (cuerpo.order ** l - 1))
        ss = Esquema(cuerpo, r, l, participantes)
    else:  # Si no es ninguno de los anteriores es el esquema simplificado
        secreto = b'\x00'
        ss = Simplificado(cuerpo, participantes)
    try:
        if orden == 'dc': # pedir decodificar antes de codificar (solo cuantico)
            ss.decodificacion([QuantumRegister(3,'a')])
        elif orden == 'ca': # pedir comparticines anticipadas despues de codificar
            ss.codificacion(secreto)
            ss.comparticion_anticipada(['a'])
        elif orden == 'cc': # pedir codificar tras haber codificado
            ss.codificacion(secreto)
            ss.codificacion(secreto)
        elif orden == 'dd': # pedir decodificar tras haber decodificado (solo cuánticos)
            part = ss.codificacion(secreto)
            ss.decodificacion(part)
            ss.decodificacion(part)
        else: # compartir anticipadamente tras compartir anticipadamente (solo cuánticos)
            ss.comparticion_anticipada(['a'])
            ss.comparticion_anticipada(['a'])
    except Exception as e:
        print(f'{type(e).__name__}: {e}')
    print()

print('#' * 80)
print(' ' * 25 + 'VALIDACION DE PARÁMETROS')
print('#' * 80)

print('=' * 80)
print('Un único participante')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_datos(Esquema, GF(2, 4), r=1, l=1, participantes=['a'])

print('=' * 80)
print('Participantes duplicados')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_datos(Esquema, GF(2, 4), r=2, l=2, participantes=['a', 'a'])

print('=' * 80)
print('Cuerpo demasiado pequeño')
print('=' * 80)
print('Esquemas de Shamir y Shamir en rampa (q = 4, n = 4)')
for Esquema in [Shamir, ShamirRampa]:
    validacion_datos(Esquema, GF(2, 2), r=3, l=2, participantes=['a', 'b', 'c', 'd'])
print('Esquema de McEliece-Sarwate (q = 7, n = 6, l = 2)')
validacion_datos(McElieceSarwate, GF(7), r=3, l=2, participantes=['a', 'b', 'c', 'd', 'e', 'f'])
print('Esquema de Cleve-Gottesman-Lo (q = 4, m = 2r-1 = 5)')
validacion_datos(CGL, GF(2, 2), r=3, l=2, participantes=['a', 'b', 'c', 'd', 'e'])
print('Esquema de Ogawa et al. (q = 4, m = 2r-l = 4)')
validacion_datos(Ogawa, GF(2, 2), r=3, l=2, participantes=['a', 'b', 'c', 'd'])
print('Esquema de Zhang-Matsumoto (q = 4, m = 2r-l = 4, l = 2)')
validacion_datos(ZhangMatsumoto, GF(2, 2), r=3, l=2, participantes=['a', 'b', 'c', 'd'])

print('=' * 80)
print('Insuficientes participantes (n = 2, r = 3)')
print('=' * 80)
for Esquema in clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_datos(Esquema, GF(2, 4), r=3, l=2, participantes=['a', 'b'])

print('=' * 80)
print('Longitud del secreto muy grande (r = 4, l = 4)')
print('=' * 80)
for Esquema in clasico_rampa + cuantico_rampa:
    validacion_datos(Esquema, GF(2, 4), r=4, l=4, participantes=['a', 'b', 'c', 'd'])

print('=' * 80)
print('Más participantes que los permitidos (n = 7, m = 2r-1 = 5 o 2r-l = 4)')
print('=' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    validacion_datos(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd', 'e', 'f', 'g'])

print('=' * 80)
print('Elemento base del cuerpo no es 2 en esquemas cuánticos')
print('=' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    validacion_datos(Esquema, GF(3, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'])

print('#' * 80)
print(' ' * 25 + 'VALIDACIÓN DEL SECRETO')
print('#' * 80)
print('=' * 80)
print('Secreto clásico demasiado grande (q = 3^40)')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral:
    validacion_secreto(Esquema, GF(3, 40), b'896tgkjp3',r=3, participantes=['a', 'b', 'c', 'd'])
for Esquema in clasico_rampa:
    validacion_secreto(Esquema, GF(3, 40), [b'748293', b'896tgkjp3'],r=3, l=2, participantes=['a', 'b', 'c', 'd'])
print('=' * 80)
print('Secreto clásico de longitud incorrecta (l = 3)')
print('=' * 80)
for Esquema in clasico_rampa:
    validacion_secreto(Esquema, GF(2, 64), [b'748293', b'896tgkj3', b'dahjka', b'hask'], r=4, l=3, participantes=['a', 'b', 'c', 'd'])
print('=' * 80)
print('Secreto cuántico de dimensión incorrecta (q = 2^3)')
print('=' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    validacion_secreto(Esquema, GF(2, 3), Statevector([1,3,1j,0,3+2j,0,6]/np.sqrt(60)), r=3, l=2, participantes=['a', 'b', 'c'])

print('=' * 80)
print('Secreto cuántico no unitario')
print('=' * 80)
for Esquema in cuantico_umbral:
    validacion_secreto(Esquema, GF(2, 2), Statevector([1, 3, 1j, 0]), r=2, participantes=['a', 'b', 'c'])
for Esquema in cuantico_rampa:
    validacion_secreto(Esquema, GF(2, 3), Statevector([-2-7j, -4+0j, -1-5j, 2j, 0, 8+7j, -9+0j, 4+5j]).tensor(Statevector([0, 8-4j, -8+5j, 0, -10-5j, 3-3j, 7+1j, -4+6j])), r=3, l=2, participantes=['a', 'b', 'c'])

print('#' * 80)
print(' ' * 25 + 'VALIDACIÓN DE PARTICIPACIONES')
print('#' * 80)
print('=' * 80)
print('Entregar participantes duplicados')
print('=' * 80)
validacion_participaciones(Simplificado, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], opcion='d')
for Esquema in clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_participaciones(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'], opcion='d')

print('=' * 80)
print('Entregar participantes inexistentes')
print('=' * 80)
validacion_participaciones(Simplificado, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], opcion='ne')
for Esquema in clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_participaciones(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'], opcion='ne')

print('=' * 80)
print('Número insuficiente de participantes (r = 4)')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_participaciones(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], opcion='i')

print('#' * 80)
print(' ' * 20 + 'VALIDACIÓN DE PARTICIPACIONES ANTICIPADAS')
print('#' * 80)
print('=' * 80)
print('Entregar participantes duplicados')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_participaciones_anticipadas(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'], opcion='d')

print('=' * 80)
print('Entregar participantes inexistentes')
print('=' * 80)
validacion_participaciones(Simplificado, GF(2, 3), 4, l=2, participantes=['a', 'b', 'c', 'd'], opcion='ne')
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_participaciones_anticipadas(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'], opcion='ne')

print('=' * 80)
print('Pedir más participaciones anticipadas de las esperadas')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_participaciones_anticipadas(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], opcion='m')

print('=' * 80)
print('Pedir participaciones para participantes a los que ya se \nha entregado participaciones anticipadas')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa:
    validacion_participaciones_anticipadas(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], opcion='df')

print('#' * 80)
print(' ' * 27 + 'VERIFICACIÓN DE ORDEN')
print('#' * 80)
print('=' * 80)
print('Pedir decodificar un secreto sin haberlo codificado anteriormente')
print('=' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    validacion_orden(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], orden='dc')

print('=' * 80)
print('Pedir compartición anticipada tras haber codificado ya el secreto')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_orden(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], orden='ca')

print('=' * 80)
print('Pedir codificar un secreto si ya se ha realizado')
print('=' * 80)
for Esquema in [Simplificado] + clasico_umbral + clasico_rampa + cuantico_umbral + cuantico_rampa:
    validacion_orden(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], orden='cc')

print('=' * 80)
print('Pedir compartición anticipada si ya se ha realizado')
print('=' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    validacion_orden(Esquema, GF(2, 3), r=4, l=2, participantes=['a', 'b', 'c', 'd'], orden='aa')

print('=' * 80)
print('Pedir decodificar si ya se ha realizado')
print('=' * 80)
for Esquema in cuantico_umbral + cuantico_rampa:
    validacion_orden(Esquema, GF(2, 3), r=3, l=2, participantes=['a', 'b', 'c', 'd'], orden='dd')