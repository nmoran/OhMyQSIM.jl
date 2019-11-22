using TensorOperations
import Base: ==, ≈
using LinearAlgebra

export QuantumRegister, FullStateQuantumRegister, to_str
export apply_1qubit_full, apply_1qubit, apply_1qubit!
export apply_2qubit_full, apply_2qubit, apply_2qubit!
export decompose_2_qubit_gate
export swap_2qubits, get_conf, measure, get_counts, measure_probs

abstract type QuantumRegister end

mutable struct FullStateQuantumRegister{T} <: QuantumRegister
    N::Integer
    state::Array{T, 1}
    FullStateQuantumRegister{T}(N) where T <: Number = new(2, zeros(T, 2^N))
    FullStateQuantumRegister{T}(N, state) where T <: Number = new(2, state)
    function FullStateQuantumRegister{T}(N, conf::String) where T <: Number
        @assert length(conf) == N
        state = 1
        for i in 1:N
            state = kron(state, conf[i] == '0' ? [1., 0] : [0, 1])
        end
        new(N, state)
    end
end

"""
    binary_repr(num, N)

Get the binary string representation
of a decimal integer.
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

function ≈(a::FullStateQuantumRegister, b::FullStateQuantumRegister)
    abs(abs(conj(a.state)' * b.state) - 1.0)  < 1e-15
end

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
    apply_1qubit(qreg, gate, i)

Apply the 1 qubit gate to qubit i inplace
"""

function apply_1qubit!(qreg::QuantumRegister, gate, i)
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
    apply_2qubit(qreg, gate, i, j)

Apply the provided 2 qubit gate to qubits i and j inplace
"""

function apply_2qubit!(qreg::QuantumRegister, gate, i, j)
    qreg.state = apply_2qubit(qreg, gate, i, j)
end

function apply_2qubit_full!(qreg::QuantumRegister, gate, i, j)
    qreg.state = apply_2qubit_full(qreg, gate, i, j)
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
