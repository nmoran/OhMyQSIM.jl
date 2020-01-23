var documenterSearchIndex = {"docs":
[{"location":"#OhMyQSIM.jl-Documentation-1","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.jl Documentation","text":"","category":"section"},{"location":"#","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.jl Documentation","text":"OhMyQSIM is a toy quantum simulator package written in Julia to experiment with different methods and approaches. If you have ideas on how to improve or extend, please get involved.","category":"page"},{"location":"#","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.jl Documentation","text":"CurrentModule = OhMyQSIM","category":"page"},{"location":"#","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.jl Documentation","text":"QuantumRegister\nFullStateQuantumRegister\nbinary_repr\nto_str\napply_1qubit_full\napply_1qubit\napply_1qubit!\ndecompose_2_qubit_gate\napply_2qubit_full\nswap_2qubits\napply_2qubit\napply_2qubit!\nmeasure_probs\nget_conf\nmeasure\nget_counts\n\n","category":"page"},{"location":"#OhMyQSIM.QuantumRegister","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.QuantumRegister","text":"Abstract type and parent of all quantum register types\n\n\n\n\n\n","category":"type"},{"location":"#OhMyQSIM.FullStateQuantumRegister","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.FullStateQuantumRegister","text":"Implementation type for full state quantum register\n\n\n\n\n\n","category":"type"},{"location":"#OhMyQSIM.binary_repr","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.binary_repr","text":"binary_repr(num, N)\n\nGet the binary string representation of a decimal integer.\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.to_str","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.to_str","text":"to_str(qreg::FullStateQuantumRegister)\n\nFunction to convert express fullstate quantum register as a string\n\nExamples\n\njulia> ψ = FullStateQuantumRegister{Int}(3, \"000\")\njulia> ψ.state[end] = 1\njulia> to_str(ψ)\n\"(1|000>) + (1|111>)\"\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.apply_1qubit_full","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.apply_1qubit_full","text":"apply_1qubit_full(qreg, gate, i)\n\nApply the 1 qubit gate to qubit i by expanding operator fully\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.apply_1qubit","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.apply_1qubit","text":"apply_1qubit(qreg, gate, i)\n\nApply the 1 qubit gate to qubit i\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.apply_1qubit!","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.apply_1qubit!","text":"apply_1qubit!(qreg, gate, i)\n\nApply the 1 qubit gate to qubit i inplace\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.decompose_2_qubit_gate","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.decompose_2_qubit_gate","text":"decompose_2_qubit_gate(gate; threshold=1e-15)\n\nDecompose a 2 qubit gate into a sum of 1 qubit gates on each qubit\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.apply_2qubit_full","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.apply_2qubit_full","text":"apply_2qubit_full(qreg, gate, i, j)\n\nApply the provided 2 qubit gate to qubits i and j by expanding full operator matrix\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.swap_2qubits","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.swap_2qubits","text":"swap_2qubits(X)\n\nGiven a two qubit operation, switch the order of the qubits\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.apply_2qubit","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.apply_2qubit","text":"apply_2qubit(qreg, gate, i, j)\n\nApply the provided 2 qubit gate to qubits i and j\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.apply_2qubit!","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.apply_2qubit!","text":"apply_2qubit!(qreg, gate, i, j)\n\nApply the provided 2 qubit gate to qubits i and j inplace\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.measure_probs","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.measure_probs","text":"measure_probs(qreg)\n\nGet for each bitstring\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.get_conf","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.get_conf","text":"get_conf(cprobs, num, N)\n\nGiven an array of cumulative probabilities of configurations and a random number from [0, 1], return the bitstring it corresponds to\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.measure","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.measure","text":"measure(qreg, shots)\n\nGiven a quantum register, return a typical set of measurements expected\n\n\n\n\n\n","category":"function"},{"location":"#OhMyQSIM.get_counts","page":"OhMyQSIM.jl Documentation","title":"OhMyQSIM.get_counts","text":"get_counts(results)\n\nGiven a set of measurement results in the form of bitstrings, return the counts of each\n\n\n\n\n\n","category":"function"}]
}