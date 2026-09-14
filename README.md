# Secret Sharing

Python implementation of several linear secret sharing schemes, including classical and quantum schemes, both threshold and ramp schemes.

In addition, all implemented schemes support advance sharing, allowing some shares to be distributed before the secret is known.

## Features

- Implementation of the Additive, Shamir, Shamir ramp and McEliece-Sarwate classical secret sharing schemes.
- Implementation of the Cleve-Gottesman-Lo, Ogawa et al., and Zhang-Matsumoto quantum secret sharing schemes.
- Built-in support for advance sharing across all schemes.
- Unified interface for all implemented algorithms.

---

## Technology Stack

The project uses the following main libraries:

| Category            | Library      |
|---------------------|--------------|
| Numerical Computing | [NumPy]      |
| Finite Fields       | [galois]     |
| Quantum Computing   | [Qiskit]     |
| Quantum Simulation  | [Qiskit Aer] |


[NumPy]:  https://github.com/numpy/numpy
[galois]: https://github.com/mhostetter/galois
[Qiskit]: https://github.com/Qiskit/qiskit
[Qiskit Aer]: https://github.com/Qiskit/qiskit-aer

The project dependencies can be installed by running:

```bash
pip install -r requirements.txt
```

---

## Using the Schemes

All scheme classes provide three main methods for performing the secret sharing procedure.

First, an instance of the corresponding scheme must be created.

```python
import galois
from classic_schemes import Shamir

scheme = Shamir(2**30, 3, ['Alice', 'Bob', 'Charles', 'Daisy'])
```

If advance sharing is desired, the `advance_sharing()` method can be called with the corresponding participants.

```python
advance = scheme.advance_sharing('Bob')
```

Regardless of whether shares have been advance shared, the `distribute()` method is used to distribute a secret.

```python
shares = scheme.distribute(b'672')
```
**Note**: If shares have been advance shared, they will **not** be returned again by the `distribute()` method.

Finally, any set of at least three participants can reconstruct the secret.

```python
secret = scheme.reconstruct([advance[0], shares[1], shares[0]])
```

### Interactive Script

An interactive script is also provided to guide the user through the execution of a scheme by requesting all required parameters through the terminal.

The script can be started with:

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