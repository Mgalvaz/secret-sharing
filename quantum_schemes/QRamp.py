import numpy as np
from galois import GF
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import partial_trace
from qiskit.circuit.library import LinearFunction
from qiskit_aer import AerSimulator

from utils import extend_matrix, simulate_statevector


class Ogawa:
    r"""
    Ogawa et al. quantum ramp secret sharing scheme over the complex Hilbert space $\mathcal{H}_{2^m}$.

    Example:
        Creates an Ogawa et al. scheme with reconstruction threshold r = 4 and secret length l = 3 over $\mathcal{H}_{2^5}$ for the participants ['a', 'b', 'c', 'd', 'e'].

        .. ipython:: python

            ogawa = Ogawa(2**5, 4, 3, ['a', 'b', 'c', 'd', 'e'])
    """
    def __init__(self, order, r, l, participants, **backend_options):
        r"""
        Creates an Ogawa et al. quantum secret sharing scheme over the complex Hilbert space $\mathcal{H}_{2^m}$.
        Quantum simulation done with AerSimulator using the statevector method.
        Other run-time options for the simulator may be specified as kwargs.
        :param order: The dimension of the Hilbert space.
        :param r: The reconstruction threshold of the scheme, i.e., the minimum number of participants required to reconstruct the secret.
        :param l: The length of the secret to be shared.
        :param participants: A list containing the unique identifiers of all participants in the scheme.
        :param backend_options: Additional options for the AerSimulator.
        """
        # Condition checks
        if order.bit_count() != 1:
            raise ValueError('The dimension of the Hilbert space must be a power of 2.')
        if 2*r-l < len(participants):
            raise ValueError(f'The number of participants ({len(participants)}) must be less than or equal to 2r - l ({2 * r - l}).')
        if order <= 2*r-l:
            raise ValueError(f'The total number of evaluation points ({2 * r - l}) must be smaller than the dimension of the Hilbert space ({order}).')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')
        if l < 2:
            raise ValueError(f'The secret length ({l}) must be greater than 1.')
        if r <= l:
            raise ValueError(f'The secret length ({l}) must be less than the reconstruction threshold ({r}).')
        if len(participants) < r:
            raise ValueError(f'The reconstruction threshold ({r}) must be less than or equal to the number of participants ({len(participants)}).')

        self.sim = AerSimulator(method='statevector', **backend_options)
        self.field = GF(order)
        self.reconstruction = r
        self.secret_length = l
        self.__advance_shares = []
        self.x_advance_remaining = []
        self.__secret_register = None
        actual_shares = [QuantumRegister(self.field.degree, name) for name in participants]
        dummy_shares = [QuantumRegister(self.field.degree, f'p{i}') for i in range(len(participants) + 1, 2 * r - l + 1)]
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
        if self.reconstruction - self.secret_length < len(advance_participants):
            raise ValueError(f'The number of advance shares ({len(advance_participants)}) must be less than or equal to the privacy threshold ({self.reconstruction - self.secret_length}).')

        qc = self.__circuit
        r = self.reconstruction
        l = self.secret_length
        x_advance = [self.participants_number[nombre] for nombre in advance_participants]
        shares_advance = [self.__all_shares[idx - 1] for idx in x_advance]
        x_not_advance = np.setdiff1d(range(1, 2 * r - l + 1), x_advance) # Numbers associated to all participants who don't receive advance shares
        shares_not_advance = [self.__all_shares[idx - 1] for idx in x_not_advance] # Registers of the participants who won't receive advance
        undistributed_advance_shares = shares_not_advance[:r - len(shares_advance) - l]
        psi1 = shares_advance + undistributed_advance_shares # Extend advance shares to have r-l
        psi2 = shares_not_advance[r - len(shares_advance) - l:-l] # Not advance shares
        self.__secret_register = shares_not_advance[-l:] # Register where the secret will be initialized
        # Create maximally entangled state (sum |y>|y>)
        for qudit_1, qudit_2 in zip(psi1, psi2):
            qc.h(qudit_1)
            qc.cx(qudit_1, qudit_2)
        self.__advance_shares = psi1
        self.x_advance_remaining = x_not_advance[:r - len(shares_advance) - l].astype(int)
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
        if secret.dim != self.field.order**self.secret_length:
            raise ValueError(f'A statevector of dimension {self.field.order ** self.secret_length} was expected, but a statevector of dimension {secret.dim} was provided.')

        qc = self.__circuit
        r = self.reconstruction
        l = self.secret_length
        # Standard procedure
        if len(self.__advance_shares) == 0:
            x = np.arange(1, 2*r-l+1)
            qc.initialize(secret, self.__all_shares[:l]) # Initialize the secret
            for share in self.__all_shares[l:r]: # Create a uniform superposition over all possible values of the polynomial coefficients
                qc.h(share)
            # Evaluate the polynomial at each participant's evaluation point
            vandermonde = self.field(x)[:, None] ** np.arange(r)
            matrix = self.field(np.column_stack([vandermonde, np.vstack([np.eye(r - l), np.zeros((r, r - l))])]))
            matrix = extend_matrix(matrix)
            participant_order = [qubit for share in self.__all_shares for qubit in reversed(share)]
        # Advance sharing
        else:
            qc.initialize(secret, self.__secret_register)
            x_advance = [self.participants_number[share.name] for share in self.__advance_shares]
            x_remaining = np.setdiff1d(range(1, 2*r-l+1), x_advance) # Numbers of all participants who don't have a share yet
            x = np.concat([x_remaining, self.x_advance_remaining])
            shares_remaining = [self.__all_shares[idx - 1] for idx in x_remaining] # Shares of all participants who don't have a share yet
            vandermonde = self.field(x_advance)[:, None] ** np.arange(r)
            matrix_s1 = np.linalg.inv(self.field(np.vstack([vandermonde, np.column_stack([np.eye(l), np.zeros((l, r - l))])])))
            matrix_s2 = self.field(x_remaining)[:, None] ** np.arange(r)
            matrix = extend_matrix(matrix_s2 @ matrix_s1)
            participant_order = [qubit for share in shares_remaining for qubit in reversed(share)]

        qc.append(LinearFunction(matrix), participant_order) # Apply evaluation matrix
        self.__advance_shares = None # Delete the stored advance shares for further security
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
        l = self.secret_length
        qc = self.__circuit
        # Obtain the evaluation points associated with the participants.
        # Sort the shares to preserve the ordering of the secret after performing the partial trace
        shares_sorted = sorted(shares[:r], key=lambda p: self.participants_number[p.name])
        x_sorted = [self.participants_number[share.name] for share in shares_sorted]
        x_remaining = np.setdiff1d(np.arange(1, 2 * r - l + 1), x_sorted)
        participant_order = [qubit for share in shares_sorted for qubit in reversed(share)]
        matrix_s1 = np.linalg.inv(self.field(x_sorted)[:, None] ** np.arange(r))
        vandermonde = self.field(x_remaining)[:,None] ** np.arange(r)
        matrix_s2 = self.field(np.vstack([np.column_stack([np.eye(l), np.zeros((l, r - l))]), vandermonde]))
        matrix = extend_matrix(matrix_s2 @ matrix_s1)
        qc.append(LinearFunction(matrix), participant_order)
        sv = simulate_statevector(qc, self.sim)
        trace_indices = np.setdiff1d(range((2*r-l) * self.field.degree), np.concatenate([range((i - 1) * self.field.degree, i * self.field.degree) for i in x_sorted[:l]])) # Position of the qubits to trace
        self.__circuit = None # Mark the reconstruction procedure as completed
        return partial_trace(sv, trace_indices.tolist()).to_statevector()

class ZhangMatsumoto:
    r"""
    Zhang-Matsumoto quantum ramp secret sharing scheme over the complex Hilbert space $\mathcal{H}_{2^m}$.

    Example:
        Creates a Zhang-Matsumoto scheme with reconstruction threshold r = 4 and secret length l = 3 over $\mathcal{H}_{2^5}$ for the participants ['a', 'b', 'c', 'd', 'e'].

        .. ipython:: python

            qss = ZhangMatsumoto(2**5, 4, 3, ['a', 'b', 'c', 'd', 'e'])
    """
    def __init__(self, order, r, l, participants, **backend_options):
        r"""
        Creates a Zhang-Matsumoto quantum ramp secret sharing scheme over the complex Hilbert space $\mathcal{H}_{2^m}$.
        :param order: The dimension of the Hilbert space.
        :param r: The reconstruction threshold of the scheme, i.e., the minimum number of participants required to reconstruct the secret.
        :param l: The length of the secret to be shared.
        :param participants: A list containing the unique identifiers of all participants in the scheme.
        :param backend_options: Additional options for the AerSimulator.
        """
        # Condition checks
        if order.bit_count() != 1:
            raise ValueError('The dimension of the Hilbert space must be a power of 2.')
        if 2*r-l < len(participants):
            raise ValueError(f'The number of participants ({len(participants)}) must be less than or equal to 2r - l ({2 * r - l}).')
        if order < 2*r:
            raise ValueError(f'The total number of evaluation points ({2 * r - l}) must be less than or equal to the order of the underlying field minus the secret length ({order - l}).')
        if len(participants) != len(set(participants)):
            raise ValueError('Duplicate participants were found.')
        if l < 2:
            raise ValueError(f'The secret length ({l}) must be greater than 1.')
        if r <= l:
            raise ValueError(f'The secret length ({l}) must be less than the reconstruction threshold ({r}).')
        if len(participants) < r:
            raise ValueError(f'The reconstruction threshold ({r}) must be less than or equal to the number of participants ({len(participants)}).')

        self.sim = AerSimulator(method='statevector', **backend_options)
        self.field = GF(order)
        self.reconstruction = r
        self.secret_length = l
        self.__advance_shares = []
        self.x_advance_remaining = []
        self.__secret_register = None
        actual_shares = [QuantumRegister(self.field.degree, name) for name in participants]
        dummy_shares = [QuantumRegister(self.field.degree, f'p{i}') for i in range(len(participants) + 1, 2 * r - l + 1)]
        self.__all_shares = actual_shares + dummy_shares
        self.__circuit = QuantumCircuit(*self.__all_shares)
        self.participants_number = {name: i for i, name in enumerate(participants, l)}

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
        if self.reconstruction - self.secret_length < len(advance_participants):
            raise ValueError(f'The number of advance shares ({len(advance_participants)}) must be less than or equal to the privacy threshold ({self.reconstruction - self.secret_length}).')

        qc = self.__circuit
        r = self.reconstruction
        l = self.secret_length
        x_advance = [self.participants_number[nombre] for nombre in advance_participants]
        shares_advance = [self.__all_shares[idx - l] for idx in x_advance]
        x_not_advance = np.setdiff1d(range(l, 2 * r), x_advance) # Numbers associated to all participants who don't receive advance shares
        shares_not_advance = [self.__all_shares[idx - l] for idx in x_not_advance] # Shares of all participants who don't receive advance shares
        undistributed_advance_shares = shares_not_advance[:r - len(shares_advance) - l]
        psi1 = shares_advance + undistributed_advance_shares # Extend advance shares to have r-l
        psi2 = shares_not_advance[r - len(shares_advance) - l:-l] # Not advance shares
        self.__secret_register = shares_not_advance[-l:] # Register where the secret will be initialized
        # Create maximally entangled state (sum |y>|y>)
        for qudit_1, qudit_2 in zip(psi1, psi2):
            qc.h(qudit_1)
            qc.cx(qudit_1, qudit_2)
        self.__advance_shares = psi1
        self.x_advance_remaining = x_not_advance[:r - len(shares_advance) - l].astype(int)
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
        if secret.dim != self.field.order**self.secret_length:
            raise ValueError(f'A statevector of dimension {self.field.order**self.secret_length} was expected, but a statevector of dimension {secret.dim} was provided.')

        qc = self.__circuit
        r = self.reconstruction
        l = self.secret_length
        # Standard procedure
        if len(self.__advance_shares) == 0:
            x = np.arange(l, 2*r)
            qc.initialize(secret, self.__all_shares[:l]) # Initialize the secret
            for share in self.__all_shares[l:r]: # Create a uniform superposition over all possible values of the remaining polynomial evaluations
                qc.h(share)
            matrix_s1 = np.linalg.inv(self.field.Range(0, r)[:, None] ** np.arange(r)) # Obtain all polynomials that evaluate A as the secret
            shares_polynomial = [qubit for share in self.__all_shares[:r] for qubit in reversed(share)]
            qc.append(LinearFunction(extend_matrix(matrix_s1)), shares_polynomial) # Apply matrix that obtains E_r(A,s)
            vandermonde = self.field(x)[:, None] ** np.arange(r) # Evaluate the polynomial at the participant evaluation points
            matrix_s2 = self.field(np.column_stack([vandermonde, np.vstack([np.eye(r - l), np.zeros((r, r - l))])])) # Evaluation matrix
            matrix = extend_matrix(matrix_s2)
            participant_order = [qubit for share in self.__all_shares for qubit in reversed(share)]
        # Advance sharing
        else:
            qc.initialize(secret, self.__secret_register)
            x_advance = [self.participants_number[share.name] for share in self.__advance_shares]
            x_remaining = np.setdiff1d(range(l, 2 * r), x_advance) # Numbers of all participants who don't have a share yet
            x = np.concat([x_remaining, self.x_advance_remaining])
            shares_remaining = [self.__all_shares[idx - l] for idx in x_remaining] # Shares of all participants who don't have a share yet
            matrix_s1 = np.linalg.inv(self.field(np.concat([x_advance, np.arange(l)]))[:, None] ** np.arange(r))
            matrix_s2 = self.field(x_remaining)[:, None] ** np.arange(r)
            matrix = extend_matrix(matrix_s2 @ matrix_s1)
            participant_order = [qubit for share in shares_remaining for qubit in reversed(share)]
        qc.append(LinearFunction(matrix), participant_order) # Apply evaluation matrix
        self.__advance_shares = None # Delete the stored advance shares for further security
        # Return the remaining actual participant shares
        return list(self.__all_shares[i - l] for i in x[x < len(self.participants_number) + l])

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
        valid_shares = set(self.__all_shares)
        for share in shares:
            if share not in valid_shares:
                raise ValueError(f"Participant '{share.name}' is not registered or provided an invalid quantum register.")

        r = self.reconstruction
        l = self.secret_length
        qc = self.__circuit
        # Sort the shares to preserve the ordering of the secret after # performing the partial trace
        shares_sorted = sorted(shares[:r], key=lambda p: self.participants_number[p.name])
        x_sorted = [self.participants_number[share.name] for share in shares_sorted]
        x_remaining = self.field(np.setdiff1d(np.arange(2 * r), x_sorted))
        participant_order = [qubit for share in shares_sorted for qubit in reversed(share)]
        matrix_s1 = np.linalg.inv(self.field(x_sorted)[:,None] ** np.arange(r))
        matrix_s2 = x_remaining[:,None] ** np.arange(r)
        matrix = extend_matrix(matrix_s2 @ matrix_s1)
        qc.append(LinearFunction(matrix), participant_order)
        sv = simulate_statevector(qc, self.sim)
        trace_indices = np.setdiff1d(range((2 * r - l) * self.field.degree), np.concatenate([range((i - l) * self.field.degree, (i - l + 1) * self.field.degree) for i in x_sorted[:l]])) # Position of the qubits to trace
        self.__circuit = None # Mark the reconstruction procedure as completed
        return partial_trace(sv, trace_indices.tolist()).to_statevector()