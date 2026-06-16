from main_clasico import programa_clasico
from main_cuantico import programa_cuantico

from utils import pedir_entero

# Pedir datos del esquema
print('¿Qué tipo de esquema desea usar?\n1.- Esquemas clásicos.\n2.- Esquemas cuánticos.')
version = pedir_entero('Respuesta: ',f'No se ha introducido un numero válido.', lambda x: 1 <= x <= 2)
if version == 1:
    programa_clasico()
else:
    programa_cuantico()