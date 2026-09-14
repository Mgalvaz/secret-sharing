import numpy as np
from galois import GF
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import partial_trace
from qiskit.circuit.library import LinearFunction
from qiskit_aer import AerSimulator

from utils import extend_matrix, simulate_statevector

class CGL:
    r"""
    Cleve-Gottesman-Lo quantum secret sharing scheme over the complex Hilbert space $\mathcal{H}_{2^m}$.

    Example:
        Creates a Cleve-Gottesman-Lo scheme with reconstruction threshold r = 4 over $\mathcal{H}_{2^5}$ for the participants ['a', 'b', 'c', 'd', 'e', 'f'].

        .. ipython:: python

            cgl = CGL(2**5, 4, ['a', 'b', 'c', 'd', 'e', 'f'])
    """
    def __init__(self, order, r, participants, **backend_options):
        r"""
        Creates a Cleve-Gottesman-Lo quantum secret sharing scheme over the complex Hilbert space $\mathcal{H}_{2^m}$.
        Quantum simulation done with AerSimulator using the statevector method.
        Other run-time options for the simulator may be specified as kwargs.
        :param order: The dimension of the Hilbert space.
        :param r: The reconstruction threshold of the scheme, i.e., the minimum number of participants required to reconstruct the secret.
        :param participants: A list containing the unique identifiers of all participants in the scheme.
        """
        # Condition checks
        if order.bit_count() != 1:
            raise ValueError('The dimension of the Hilbert space must be a power of 2.')
        if 2*r-1 < len(participants):
            raise ValueError(f'The number of participants ({len(participants)}) must be less than or equal to 2r - 1 ({2 * r - 1}).')
        if order <= 2*r-1:
            raise ValueError(f'The total number of evaluation points ({2 * r - 1}) must be smaller than the dimension of the Hilbert space ({order}).')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')
        if r < 2:
            raise ValueError(f'The reconstruction threshold ({r}) must be greater than 1.')
        if len(participants) < r:
            raise ValueError(f'The reconstruction threshold ({r}) must be less than or equal to the number of participants ({len(participants)}).')

        self.sim = AerSimulator(method='statevector', **backend_options)
        self.field = GF(order)
        self.reconstruction = r
        self.__advance_shares = []
        self.x_advance_remaining = []
        self.__secret_register = None
        actual_shares = [QuantumRegister(self.field.degree, participant) for participant in participants]
        dummy_shares = [QuantumRegister(self.field.degree, f'p{i}') for i in range(len(participants) + 1, 2 * r)]
        self.__all_shares = actual_shares + dummy_shares
        self.__circuit = QuantumCircuit(*self.__all_shares)
        self.participants_number = {name: i for i, name in enumerate(participants, 1)}

    def advance_sharing(self, advance_participants):
        """
        Creates advance shares for the specified participants.
        Each share is represented by a quantum register whose name corresponds to the unique identifier of the participant receiving it.
        Due to quantum limitations, this method can only be called once.
        :param advance_participants: Sequence of participants that will receive advance shares.
        :return: A list containing the quantum registers assigned to the specified participants.
        """
        # Condition checks
        if self.__advance_shares is None:
            raise AttributeError('All shares have already been distributed.')
        if self.__secret_register is not None:
            raise AttributeError('Advance sharing has already been performed.')
        if len(advance_participants) != len(set(advance_participants)):
            raise ValueError('Duplicate participants were found.')
        for name in advance_participants:
            if name not in self.participants_number:
                raise ValueError(f"Participant '{name}' is not registered.")
        if self.reconstruction <= len(advance_participants):
            raise ValueError(f'The number of advance shares ({len(advance_participants)}) must be less than or equal to the privacy threshold ({self.reconstruction - 1}).')

        qc = self.__circuit
        r = self.reconstruction
        x_advance = [self.participants_number[nombre] for nombre in advance_participants]
        shares_advance = [self.__all_shares[idx - 1] for idx in x_advance]
        x_not_advance = np.setdiff1d(range(1, 2 * r), x_advance) # Numbers associated to all participants who don't receive advance shares
        shares_not_advance = [self.__all_shares[idx - 1] for idx in x_not_advance] # Shares of all participants who don't receive advance shares
        undistributed_advance_shares = shares_not_advance[:r - len(shares_advance) - 1]
        psi1 = shares_advance + undistributed_advance_shares # Extend advance shares to have r-1
        psi2 = shares_not_advance[r - len(shares_advance) - 1:-1] # Not advance shares
        self.__secret_register = shares_not_advance[-1] # Register where the secret will be initialized
        # Create maximally entangled state (sum |y>|y>)
        for qudit_1, qudit_2 in zip(psi1, psi2):
            qc.h(qudit_1)
            qc.cx(qudit_1, qudit_2)
        self.__advance_shares = psi1
        self.x_advance_remaining = x_not_advance[:r - len(shares_advance) - 1].astype(int)
        return shares_advance

    def distribute(self, secret):
        """
        Creates the shares for all participants according to the given secret.
        Each share is represented by a quantum register whose name corresponds to the unique identifier of the participant receiving it.
        If advance sharing has been done, the generated shares will be consistent with advance shares.
        :param secret: The quantum secret to be shared among the participants.
        :return: A list containing the shares of all participants that did not receive an advance share.
        """
        # Condition checks
        if self.__advance_shares is None:
            raise AttributeError('All shares have already been distributed.')
        if not secret.is_valid():
            raise ValueError('A valid quantum state was not provided.')
        if secret.dim != self.field.order:
            raise ValueError(f'A statevector of dimension {self.field.order} was expected, but a statevector of dimension {secret.dim} was provided.')

        qc = self.__circuit
        r = self.reconstruction
        # Standard procedure
        if len(self.__advance_shares) == 0:
            x = np.arange(1, 2*r)
            qc.initialize(secret, [self.__all_shares[0]])  # Initialize the secret
            for share in self.__all_shares[1:r]:  # Create a uniform superposition over all possible values of the polynomial coefficients
                qc.h(share)
            # Evaluate the polynomial at each participant's evaluation point
            vandermonde = self.field(x)[:, None] ** np.arange(r)
            matrix = self.field(np.column_stack([vandermonde, np.vstack([np.eye(r - 1), np.zeros((r, r - 1))])])) # Evaluation matrix
            matrix = extend_matrix(matrix)  # Extend matrix from F_q numbers to F_2 vectors
            participant_order = [qubit for share in self.__all_shares for qubit in reversed(share)]  # Qiskit uses little-endian ordering, whereas the extended matrix is represented in big-endian order.
        # Advance sharing
        else:
            qc.initialize(secret, self.__secret_register)
            x_advance = [self.participants_number[share.name] for share in self.__advance_shares]
            x_remaining = np.setdiff1d(range(1, 2*r), x_advance) # Numbers of all participants who don't have a share yet
            x = np.concat([x_remaining, self.x_advance_remaining])
            x_advance = self.field(x_advance + [0])
            shares_remaining = [self.__all_shares[idx - 1] for idx in x_remaining] # Shares of all participants who don't have a share yet
            matrix_s1 = np.linalg.inv(x_advance[:, None] ** np.arange(r))
            matrix_s2 = self.field(x_remaining)[:, None] ** np.arange(r)
            matrix = extend_matrix(matrix_s2 @ matrix_s1)
            participant_order = [qubit for share in shares_remaining for qubit in reversed(share)]

        qc.append(LinearFunction(matrix), participant_order)  # Apply evaluation matrix
        self.__advance_shares = None  # Delete the stored advance shares for further security
        # Return the remaining actual participant shares
        return list(self.__all_shares[i - 1] for i in x[x <= len(self.participants_number)])

    def reconstruct(self, shares):
        """
        Reconstructs the secret encoded in the provided shares.
        Each share is represented by a quantum register whose name corresponds to the unique identifier of the participant providing it.
        :param shares: Sequence containing the shares provided by the participants wishing to reconstruct the secret.
        :return: The reconstructed secret up to a global phase.
        """
        # Condition checks
        if self.__circuit is None:
            raise AttributeError('The reconstruction procedure has already been performed.')
        if self.__advance_shares is not None:
            raise AttributeError('The sharing procedure has not yet been performed.')
        if len(shares) < self.reconstruction:
            raise ValueError('Not enough shares were provided to recover the secret.')
        if len(shares) != len(set(shares)):
            raise ValueError('Duplicate participants were found.')
        shares_set = set(self.__all_shares)
        for share in shares:
            if share not in shares_set:
                raise ValueError(f"Participant '{share.name}' is not registered or provided an invalid quantum register.")

        r = self.reconstruction
        qc = self.__circuit
        # Obtain the evaluation points associated with the participants
        x_reconstruct = self.field([self.participants_number[share.name] for share in shares[:r]])
        x_remaining = self.field(np.setdiff1d(np.arange(2 * r), x_reconstruct))
        participant_order = [qubit for share in shares[:r] for qubit in reversed(share)]
        matrix_s1 = np.linalg.inv(x_reconstruct[:,None] ** np.arange(r)) # Reconstruction procedure first step matrix
        matrix_s2 = x_remaining[:,None] ** np.arange(r) # Reconstruction procedure second step matrix
        matrix = extend_matrix(matrix_s2 @ matrix_s1)
        qc.append(LinearFunction(matrix), participant_order) # Perform both steps at the same time
        sv = simulate_statevector(qc, self.sim)
        trace_indices = list(range((int(x_reconstruct[0]) - 1) * self.field.degree)) + list(range(int(x_reconstruct[0]) * self.field.degree, (2 * r - 1) * self.field.degree))  # Position of the qubits to trace
        self.__circuit = None  # Mark the reconstruction procedure as completed
        return partial_trace(sv, trace_indices).to_statevector()