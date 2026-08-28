from classical_main import programa_clasico
from quantum_main import programa_cuantico

from utils import ask_int

# Pedir datos del esquema
print('¿Qué tipo de esquema desea realizar?\n1.- Esquemas clásicos.\n2.- Esquemas cuánticos.')
version = ask_int('Respuesta: ', f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
if version == 1:
    programa_clasico()
else:
    programa_cuantico()