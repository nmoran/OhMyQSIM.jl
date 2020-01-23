export qft_circuit

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
