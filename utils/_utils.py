import secrets
import numpy as np
from base64 import b64decode, b64encode
from galois import Poly

def bytes_to_int(string):
    """
    Returns the integer representation of a bytes object.
    If a sequence of byte strings is provided, returns a list containing the integer representation of each byte string.
    :param string: The byte string whose integer representation is required.
    :return: The integer/s representation of the provided byte string/sequence.
    """
    if isinstance(string, bytes):
        return int.from_bytes(string, byteorder='big')
    return list(int.from_bytes(b, byteorder='big') for b in string)

def int_to_bytes(num):
    """
    Returns the byte representation of an integer.
    If a sequence of integers is provided, returns a list containing the byte representation of each integer.
    :param num: The integer whose byte representation is required.
    :return: The byte representation of the provided integer.
    """
    if isinstance(num, int):
        length = (num.bit_length() + 7) // 8 if num > 0 else 1
        return num.to_bytes(length, byteorder='big')
    elif isinstance(num, np.ndarray):
        if num.ndim == 0:
            n_int = int(num)
            length = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            return n_int.to_bytes(length, byteorder='big')
        representations = []
        for n_int in num.tolist():
            length = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            representations.append(n_int.to_bytes(length, byteorder='big'))
        return representations
    else:
        representations = []
        for n_int in map(int, num):
            length = (n_int.bit_length() + 7) // 8 if n_int > 0 else 1
            representations.append(n_int.to_bytes(length, byteorder='big'))
        return representations

def int_to_b64str(array, length):
    """
    Returns the Base64 encoding of each integer in an integer array.
    :param array: Integer array.
    :param length: The number of bytes used to represent each integer. An OverflowError is raised if an integer cannot be represented using the specified number of bytes.
    :return: A list containing the Base64 encoding of each integer.
    """
    if isinstance(array, np.ndarray):
        return list(b64encode(num.to_bytes(length, byteorder='big')).decode() for num in array.tolist())
    return list(b64encode(num.to_bytes(length, byteorder='big')).decode() for num in array)

def b64str_to_int(array):
    """
    Returns the integer represented by each Base64-encoded string.
    :param array: The list of encoded strings.
    :return: A list containing the decoded integers.
    """
    return list(int.from_bytes(b64decode(str_b64), byteorder='big') for str_b64 in array)

def random_array(sup, n):
    """
    Creates an array of length n containing cryptographically secure random integers in the range [0, sup).
    :param sup: Upper bound (exclusive) for the random integers.
    :param n: The number of random integers to generate.
    :return: The array of random integers.
    """
    return [secrets.randbelow(sup) for _ in range(n)]

def random_polynomial(field, degree):
    """
    Creates a cryptographically secure random polynomial over a field.
    :param field: The field over which the polynomial is constructed.
    :param degree: The degree of the polynomial. Equivalently, the number of polynomial coefficients minus one.
    :return: The random polynomial.
    """
    return Poly([secrets.randbelow(field.order) for _ in range(degree+1)], field=field)

def extend_matrix(matrix):
    """
    Returns the extension of the given matrix obtained by replacing each of its elements with its matrix representation over the base field and converting the result to Boolean format.
    :param matrix: A matrix over the field GF(2, m) to be extended.
    :return: The extended matrix over the base field GF(2).
    """
    gf = type(matrix)  # Get the working field
    if gf.characteristic != 2:
        raise ValueError(f'A field with characteristic 2 was expected, but {gf.characteristic} was provided.')
    num_bits = gf.degree
    final_shape = np.array(matrix.shape) * num_bits # Obtain the final shape
    basis = gf.primitive_element ** np.arange(num_bits - 1, -1, -1) # Basis of the extended field, represented in Big Endian order
    new_basis_matrix = gf(np.kron(matrix, basis)) # Convert the matrix to the new basis
    matriz_gf2 = new_basis_matrix.vector() # Obtain the representation of each matrix element in the new field
    matriz_gf2 = np.transpose(matriz_gf2, (0, 2, 1)).reshape(final_shape) # Align dimensions
    return matriz_gf2.view(np.ndarray).astype(bool) # Return the extended matrix in Boolean format

def ask_int(question, error_message, condition):
    """
    Prompts the user to enter an integer and validates its value.
    :param question: Message displayed to prompt the user for the integer.
    :param error_message: Message displayed if the integer does not satisfy the condition.
    :param condition: Function that receives the entered integer and returns True if it is valid and False otherwise.
    :return: The validated integer entered by the user.
    """
    while True:
        try:
            num = int(input(question))
        except ValueError:
            print('A valid number was not entered.')
        else:  # Executed if no exception occurred
            if not condition(num):
                print(error_message)
            else:
                print() # Print a blank line for better visual clarity
                return num