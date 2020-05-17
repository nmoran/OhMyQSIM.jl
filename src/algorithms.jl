export qft_circuit, simple_state_prep_circuit, ghz_circuit

"""
    qft_circuit(n::Integer)

Creates a QFT circuit
"""
function qft_circuit(n::Integer)
    qc = Gates.QuantumCircuit(n)
    for i in 1:n
        push!(qc.ops, Gates.GateOperator([i], Gates.h))
        for j in i:n-1
            push!(qc.ops, Gates.GateOperator([j+1 , i], Gates.rk(j-i+2)))
        end
    end
    for i = 1:convert(Integer, floor(n//2))
        push!(qc.ops, Gates.GateOperator([i , n + 1 - i], Gates.swap))
    end
    qc
end

"""
    function simple_state_prep_circuit(qubits::Integer, depth::Integer)

Create a simple state preparation circuit with depth layers
"""
function simple_state_prep_circuit(qubits::Integer, depth::Integer)
    qc = Gates.QuantumCircuit(qubits)
    for d in 1:depth
        angles = rand(qubits, 3) * 2π
        for i in 1:qubits
            gate = Gates.u3(angles[i, :]...)
            push!(qc.ops, Gates.GateOperator([i], gate))
        end
        for i in ((d % 2)+1):2:qubits-1
            push!(qc.ops, Gates.GateOperator([i , i+1], Gates.cx))
        end
    end
    qc
end

"""
    ghz_circuit(n::Integer)

Creates a GHZ preparation circuit
"""
function ghz_circuit(n::Integer)
    qc = Gates.QuantumCircuit(n)
    push!(qc.ops, Gates.GateOperator([1], Gates.h))
    for i in 1:n-1
        push!(qc.ops, Gates.GateOperator([i, i+1], Gates.cx))
    end
    qc
end
