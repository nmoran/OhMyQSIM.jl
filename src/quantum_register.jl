using TensorOperations
import Base: ==, ≈
using LinearAlgebra

export QuantumRegister, FullStateQuantumRegister, to_str
export apply_1qubit_full, apply_1qubit, apply_1qubit!
export apply_2qubit_full, apply_2qubit, apply_2qubit!
export decompose_2_qubit_gate
export swap_2qubits, get_conf, measure, get_counts, measure_probs
export binary_repr, execute!

"Abstract type and parent of all quantum register types"
abstract type QuantumRegister end

"Implementation type for full state quantum register"
mutable struct FullStateQuantumRegister{T} <: QuantumRegister
    N::Integer
    state::Array{T, 1}
    FullStateQuantumRegister{T}(N) where T <: Number = new(NTuple, zeros(T, 2^N))
    FullStateQuantumRegister{T}(N, state) where T <: Number = new(N, state)
    function FullStateQuantumRegister{T}(N, conf::String) where T <: Number
        @assert length(conf) == N
        new(N, configuration2vector(conf))
    end
end

FullStateQuantumRegister(N) = FullStateQuantumRegister{ComplexF64}(N)
FullStateQuantumRegister(N, state) = FullStateQuantumRegister{ComplexF64}(N, state)
FullStateQuantumRegister(N, conf::String) = FullStateQuantumRegister{ComplexF64}(N, conf)

"""
    configuration2vector(conf::String)

Convert a configuration given as a string to a hilbert space vector
# Examples
```jldoctest
julia> configuration2vector("01")
4-element Array{Float64,1}:
 0.0
 1.0
 0.0
 0.0
```
"""
function configuration2vector(conf::String)
    state = 1
    for i in 1:length(conf)
        state = kron(state, conf[i] == '0' ? [1., 0.] : [0., 1.])
    end
    state
end

"""
    binary_repr(num, N)

Get the binary string representation
of a decimal integer.
# Examples
```jldoctest
julia> println(binary_repr(5, 6))
"000101"
```
"""
function binary_repr(num, N)
    if num >= 1 << N
        num %= 1 << N
    end
    binary_str = ""
    while N > 0
        if num >= (1 << (N - 1))
            binary_str *= "1"
            num -= 1 << (N - 1)
        else
            binary_str *= "0"
        end
        N -= 1
    end
    return binary_str
end

"""
    to_str(qreg::FullStateQuantumRegister)

Function to convert express fullstate quantum register as a string

# Examples
```jldoctest
julia> ψ = FullStateQuantumRegister{Int}(3, "000")
julia> ψ.state[end] = 1
julia> to_str(ψ)
"(1|000>) + (1|111>)"
```
"""
function to_str(qreg::FullStateQuantumRegister)
    N = qreg.N
    str_rep = ""
    for i in 1:length(qreg.state)
        if abs(qreg.state[i]) > 1e-15
            if str_rep != ""
                str_rep *= " + "
            end
            str_rep *= "($(qreg.state[i])|$(binary_repr(i-1, N))>)"
        end
    end
    str_rep
end

"""
    ≈(a::FullStateQuantumRegister{T},
      b::FullStateQuantumRegister{T}) where T <: Complex

Function which does an approximate comparision of quantum registers with complex elements
"""
function ≈(a::FullStateQuantumRegister{T},
           b::FullStateQuantumRegister{T}) where T <: Complex
    maximum(abs.(a.state - b.state)) < eps(T.types[1])*log2(length(a.state))*2
end

"""
    ≈(a::FullStateQuantumRegister{T},
      b::FullStateQuantumRegister{T}) where T <: Complex

Function which does an approximate comparision of quantum registers with float elements
"""
function ≈(a::FullStateQuantumRegister{T},
           b::FullStateQuantumRegister{T}) where T <: AbstractFloat
    maximum(abs.(a.state - b.state)) < eps(T)*log2(length(a.state))*2
end

"""
    ==(a::FullStateQuantumRegister{T},
       b::FullStateQuantumRegister{T})

Function which does an exact comparision of elements of quantum registers
"""
function ==(a::FullStateQuantumRegister, b::FullStateQuantumRegister)
    a.state == b.state
end

"""
    apply_1qubit_full(qreg, gate, i)

Apply the 1 qubit gate to qubit i by expanding operator fully
"""
function apply_1qubit_full(qreg::QuantumRegister, gate, i)
    op = 1
    for j = 1:qreg.N
        op = kron(op, j == i ? gate : Gates.Id)
    end
    op * qreg.state
end

"""
    apply_1qubit(qreg, gate, i)

Apply the 1 qubit gate to qubit i
"""
function apply_1qubit(qreg::QuantumRegister, gate, i)
    # reshape to separate link to apply tensor to
    state = reshape(qreg.state, (2^(qreg.N-i), 2, 2^(i-1)))
    @tensor state[a, d, c] := gate[d, b] * state[a, b, c]
    reshape(state, (2^qreg.N))
end

"""
    apply_1qubit!(qreg, gate, i)

Apply the 1 qubit gate to qubit i inplace
"""
function apply_1qubit!(qreg::QuantumRegister, gate, i)
    # reshape to separate link to apply tensor to
    # state = reshape(qreg.state, (2^(qreg.N-i), 2, 2^(i-1)))
    #@tensor state[a, d, c] = gate[d, b] * state[a, b, c]
    # state = reshape(state, (2^qreg.N))
    qreg.state = apply_1qubit(qreg, gate, i)
end


"""
    decompose_2_qubit_gate(gate; threshold=1e-15)

Decompose a 2 qubit gate into a sum of 1 qubit gates on each qubit
"""
function decompose_2_qubit_gate(gate; threshold=1e-15)
    gate = reshape(gate, (2, 2, 2, 2))
    gate = permutedims(gate, (1 ,3, 2, 4))
    gate = reshape(gate, (4, 4))
    F = svd(gate)
    s = [x for x in F.S if x > threshold]
    ops = Array{Any, 1}()
    for i in 1:length(s)
        push!(ops, (reshape(F.U[:, i:i] * s[i], (2, 2)), reshape(F.Vt[i:i, :], (2, 2))))
    end
    ops
end


"""
    apply_2qubit_full(qreg, gate, i, j)

Apply the provided 2 qubit gate to qubits i and j by expanding full operator matrix
"""
function apply_2qubit_full(qreg::QuantumRegister, gate, i, j)
    gate_ops = decompose_2_qubit_gate(gate)
    total_op = nothing
    for gate_op in gate_ops
        op = 1
        for k = 1:qreg.N
            op = kron(op, k == j ? gate_op[1] : (k == i ? gate_op[2] : Gates.Id))
        end
        if total_op == nothing
            total_op = op
        else
            total_op += op
        end
    end
    total_op * qreg.state
end

"""
    apply_2qubit_full!(qreg, gate, i, j)

Apply the provided 2 qubit gate to qubits i and j by expanding full operator matrix
"""
function apply_2qubit_full!(qreg::QuantumRegister, gate, i, j)
    qreg.state = apply_2qubit_full(qreg, gate, i, j)
end

"""
    swap_2qubits(X)

Given a two qubit operation, switch the order of the qubits
"""
function swap_2qubits(X)
    X = reshape(X, (2, 2, 2, 2))
    X = permutedims(X, (2, 1, 4, 3))
    reshape(X, (4, 4))
end

"""
    apply_2qubit(qreg, gate, i, j)

Apply the provided 2 qubit gate to qubits i and j
"""
function apply_2qubit(qreg::QuantumRegister, gate, i, j)
    if i > j
        return apply_2qubit(qreg, swap_2qubits(gate), j, i)
    else
        # reshape to separate i and j indices
        state = reshape(qreg.state, (2^(qreg.N-j), 2, 2^(j-i-1), 2, 2^(i-1)))
        gate = reshape(gate, (2, 2, 2, 2))
        @tensor state[a, b, c, d, e] := gate[b, d, x, y] * state[a, x, c, y, e]
        reshape(state, (2^qreg.N))
    end
end

"""
    apply_2qubit!(qreg, gate, i, j)

Apply the provided 2 qubit gate to qubits i and j inplace
"""
function apply_2qubit!(qreg::QuantumRegister, gate, i, j)
    qreg.state = apply_2qubit(qreg, gate, i, j)
end

"""
    measure_probs(qreg)

Get for each bitstring
"""
function measure_probs(qreg::QuantumRegister)
    real(qreg.state .* conj(qreg.state))
end

"""
    get_conf(cprobs, num, N)

Given an array of cumulative probabilities of
configurations and a random number from [0, 1],
return the bitstring it corresponds to
"""
function get_conf(cprobs, num, N)
    binary_repr(searchsortedfirst(cprobs, num)-1, N)
end

"""
    measure(qreg, shots)

Given a quantum register, return a typical set of
measurements expected
"""
function measure(qreg::QuantumRegister; shots=100)
    cprobs = cumsum(measure_probs(qreg))
    rands = rand(shots)
    return [get_conf(cprobs, x, qreg.N) for x in rands]
end

"""
    get_counts(results)

Given a set of measurement results in the form of bitstrings,
return the counts of each
"""
function get_counts(results)
    counts = Dict{String, Integer}()
    for conf in results
        if haskey(counts, conf)
            counts[conf] += 1
        else
            counts[conf] = 1
        end
    end
    counts
end

"""
    execute!(qc::QuantumCircuit, state::FullStateQuantumRegister)

Apply the given quantum circuit to the given quantum register in place
"""
function execute!(qc::Gates.QuantumCircuit, state::FullStateQuantumRegister)
    for op in qc.ops
        if op.n == 1
            apply_1qubit!(state, op.gate, op.qubits[1])
        elseif op.n == 2
            apply_2qubit!(state, op.gate, op.qubits[1], op.qubits[2])
        end
    end
end
