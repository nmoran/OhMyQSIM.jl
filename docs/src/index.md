# Getting Started

OhMyQSIM is a toy quantum simulator package written in Julia to experiment with different methods and approaches. 
If you have ideas on how to improve or extend, please get involved. A good place to start is by looking through the issues
and adding or trying to address them as appropriate.

## Notation

Quantum registers with n-qubits have qubits labeled from 1 to n from left to right as

``| m_1 m_2 m_3 \rangle = |m_1\rangle \otimes |m_2\rangle \otimes |m_3\rangle``

This mens the index of such a state in a state vector is
``m_1 2^2 + m_2 2^1 + m_3 2^0``


## Basic usage

```
using OhMyQSIM

# Create a quantum register with 3 qubits and initialise to |000>
qreg = FullStateQuantumRegister{ComplexF64}(3, "000")

# Apply a Pauli x operation to the first qubit 
qreg.apply_1_qubit!(qreg, Gates.x, 1)

# Print the register
println(to_str(qreg))
```

results in

```
(1.0 + 0.0im|100>)
```

## Basic functions

```@meta
CurrentModule = OhMyQSIM
```

```@docs
QuantumRegister
FullStateQuantumRegister
binary_repr
to_str
apply_1qubit_full
apply_1qubit
apply_1qubit!
decompose_2_qubit_gate
apply_2qubit_full
swap_2qubits
apply_2qubit
apply_2qubit!
measure_probs
get_conf
measure
get_counts


```
