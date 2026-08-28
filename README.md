# Compartición de secretos

Implementación en Python de diversos esquemas lineales de compartición de secretos, incluyendo esquemas clásicos y cuánticos, tanto de umbral como en rampa.

Además, todos los esquemas implementados incluyen la variante de compartición anticipada, que permite distribuir algunas participaciones antes de conocer el secreto.

## Características

- Implementación de los esquemas clásicos de Shamir, Shamir en rampa, McEliece-Sarwate y una versión simplificada para los esquemas de umbral.
- Implementación de los esquemas cuánticos de Cleve-Gottesman-Lo, Ogawa et al. y Zhang-Matsumoto.
- Soporte para compartición anticipada integrado en todos los esquemas.
- Interfaz unificada para todos los algoritmos implementados.

---

## Stack tecnológico

El proyecto utiliza las siguientes librerías principales:

| Categoría           | Librería     |
|---------------------|--------------|
| Cálculo numérico    | [NumPy]      |
| Cuerpos finitos     | [galois]     |
| Circuitos Cuánticos | [Qiskit]     |
| Simulación cuántica | [Qiskit Aer] |


[NumPy]:  https://github.com/numpy/numpy
[galois]: https://github.com/mhostetter/galois
[Qiskit]: https://github.com/Qiskit/qiskit
[Qiskit Aer]: https://github.com/Qiskit/qiskit-aer

La instalación de las dependencias se puede realizar ejecutando el siguiente comando.

```bash
pip install -r requirements.txt
```

---

## Realización de los esquemas

Todas las clases creadas cuentan con tres métodos para realizar los esquemas. 

Antes de nada se debe construir un objeto de la clase correspondiente.

```python
import galois
from classic_schemes import Shamir

scheme = Shamir(galois.GF(2, 30), 3, ['Alice', 'Bob', 'Charles', 'Daisy'])
```

Si se desea realizar compartición anticipada, se debe llamar al método `comparticion_anticipada` con los participantes correspondientes.
```python
advance = scheme.comparticion_anticipada('Bob', 'Daisy')
```

Independientemente de si se ha realizado la compartición anticipada o no, en caso de querer compartir un secreto se debe llamar al método `codificacion`.
```python
shares = scheme.codificacion(b'672')
```
**Nota:** En caso de haber distribuido participaciones de forma anticipada, estas **no** volverán a devolverse durante la codificación.

Finalmente, cualquier conjunto de tres o más participantes puede reconstruir el secreto.
```python
secreto = scheme.decodificacion([advance[0], shares[1], shares[0]])
```

### _Script_ interactivo

También se proporciona un _script_ interactivo que guía al usuario durante la ejecución del esquema, solicitando por terminal todos los parámetros necesarios.

Este _script_ se puede iniciar mediante el siguiente comando.

```bash
python main.py
```

---

## Project structure
```
secret-sharing
├── classical_schemes
│   ├── __init__.py
│   ├── Perfect.py
│   └── Ramp.py
├── quantum_schemes
│   ├── __init__.py
│   ├── QPerfect.py
│   └── QRamp.py
├── utils
│   ├── __init__.py
│   └── _utils.py
├── classical_main.py
├── main.py
├── quantum_main.py
├── README.md
└── requirements.txt
```