import Base.convert
export MPSQuantumRegister, convert, two_qubit_gate_to_mpo, mpo_to_two_qubit_gate
export print_info


"Implementation type MPS quantum register"
mutable struct MPSQuantumRegister{T} <: QuantumRegister
    N::Integer
    state::Array{Array{T, 3}, 1}
    function MPSQuantumRegister{T}(N::Integer, conf::String) where T <: Number
        @assert length(conf) == N
        new(N, [init_mps_tensor(c, T) for c in conf])
    end
end

MPSQuantumRegister(N::Integer, conf::String) = MPSQuantumRegister{ComplexF64}(N, conf)

"""
    convert(::Type{FullStateQuantumRegister}, x::MPSQuantumRegister)

Convert an MPS quantum register to a full state quantum register
"""
function convert(::Type{FullStateQuantumRegister{T}}, x::MPSQuantumRegister{U}) where T <: Number where U <: Number
    state = ones(U, 1, 1)
    for i in x.N:-1:1
        @tensor state[a, b, c] := state[a, j] * x.state[i][b, c, j]
        dims = size(state)
        state = reshape(state, (dims[1] * dims[2], dims[3]))
    end
    dims = size(state)
    state = reshape(state, (dims[1] * dims[2]))
    return FullStateQuantumRegister{T}(x.N, state)
end

"""
    init_mps_tensor(qubit_value::Char, ::Type{T})

Initialise an array of MPS tensors from a string configuration
"""
function init_mps_tensor(qubit_value::Char, ::Type{T}) where T <: Number
    if qubit_value == '0'
        state = [1., 0.]
    else
        state = [0, 1]
    end
    return reshape(convert(Array{T}, state), (2, 1, 1)) # idx, up, down
end

function to_str(qreg::MPSQuantumRegister{T}) where T <: Number
    to_str(convert(FullStateQuantumRegister{T}, qreg))
end

"""
    apply_1qubit!(qreg::MPSQuantumRegister, gate, i)

Apply the 1 qubit gate to qubit i
"""
function apply_1qubit!(qreg::MPSQuantumRegister, gate, i)
    @tensor qreg.state[i][a, b, c] := qreg.state[i][m, b, c] * gate[a, m]
end

"""
    two_qubit_gate_to_mpo(gate::Array{T, 2}, n::Integer) where T <: Number

Convert a gate into an mpo
"""
function two_qubit_gate_to_mpo(gate::Array{T, 2}, threshold::AbstractFloat=1e-15) where T <: Number
    gate = reshape(gate, (2, 2, 2, 2))  #  r2, r1, l2, l1
    gate = permutedims(gate, (2, 4, 1, 3)) # r1, l1, r2, l2
    gate = reshape(gate, (4, 4)) # (r1, l1), (r2, l2)
    F = svd(gate)
    mpo = Array{Array{T, 4}, 1}()
    chi = sum(F.S .> threshold)
    push!(mpo, reshape(F.U[:,1:chi] * Diagonal(F.S[1:chi]), (2, 2, 1, chi))) # r1, l1, Up, Down
    push!(mpo, permutedims(reshape(F.Vt[1:chi,:], (chi, 2, 2, 1)), (2, 3, 1, 4))) # r2, l2, Up, Down
    mpo
end

"""
    mpo_to_two_qubit_gate(mpo::Array{Array{T, 4}, 1}) where T <: Number

Convert a gate into an mpo
"""
function mpo_to_two_qubit_gate(mpo::Array{Array{T, 4}, 1}) where T <: Number
    a = reshape_tensor(mpo[1], (1, 2, (3, 4))) # r1, l1, down
    b = reshape_tensor(mpo[2], (1, 2, (3, 4))) # r2, l2, up
    @tensor c[a, b, c, d] := a[a, b, x] * b[c, d, x] # r1, l1, r2, l2
    reshape(permutedims(c, (3, 1, 4, 2)), (4, 4))
end

"""
    apply_2qubit!(qreg::MPSQuantumRegister, gate, i, j)

Apply the 2 qubit gate to qubits i and j, assumes |i - j| = 1
"""
function apply_2qubit!(qreg::MPSQuantumRegister, gate, i, j)
    @assert abs(i -j) == 1
    if j < i
        apply_2qubit!(qreg, swap_2qubits(gate), j, i)
    else
        # reshape gate to mpo
        gate_mpo = two_qubit_gate_to_mpo(gate)
        # remove upward facing dim 1 index to simplify contraction
        gate_mpo_1 = reshape_tensor(gate_mpo[1], [[1, 3], 2, 4])

        # remove downward facing dim 1 index to simplify contraction
        gate_mpo_2 = reshape_tensor(gate_mpo[2], [[1, 4], 2, 3])

        # tensor gate mpo with state mpos
        @tensor state_1[qubit_i, up1, down1, down2] :=  qreg.state[i][qubit_i_c, up1, down1] * gate_mpo_1[qubit_i, qubit_i_c, down2]
        @tensor state_2[qubit_j, up1, up2, down1] :=  qreg.state[j][qubit_j_c, up1, down1] * gate_mpo_2[qubit_j, qubit_j_c, up2]

        # reshape
        qreg.state[i] = reshape_tensor(state_1, [1, 2, [3, 4]])
        qreg.state[j] = reshape_tensor(state_2, [1, [2, 3], 4])
    end
    return
end

function print_info(qreg::MPSQuantumRegister)
    for i = 1:qreg.N
        dims = size(qreg.state[i])[2:3]
        println("$(i) has virtual dims $(dims)")
    end
end
