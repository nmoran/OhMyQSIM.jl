import Base.convert
export MPSQuantumRegister, convert, two_qubit_gate_to_mpo, mpo_to_two_qubit_gate
export print_info, compress!, enaglemment_entropy, execute!


"Implementation type MPS quantum register"
struct MPSQuantumRegister{T} <: QuantumRegister
    N::Integer
    state::Array{Array{T, 3}, 1}
    s_values::Array{Union{Nothing, Array{<:AbstractFloat, 1}}, 1}
    function MPSQuantumRegister{T}(N::Integer, conf::String) where T <: Number
        @assert length(conf) == N
        new(N, [init_mps_tensor(c, T) for c in conf], fill(nothing, N-1))
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
    # @assert abs(i -j) == 1
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

        if abs(i -j) > 1
            chi = size(gate_mpo_1)[3]
            filler_mpo = tensor_identity(chi, chi, 2)
            for k in i+1:j-1
                @tensor tmp[qubit_k, up1, up2, down1, down2] :=
                            qreg.state[k][qubit_k_c, up1, down1] *
                            filler_mpo[up2, down2, qubit_k_c, qubit_k]
                qreg.state[k] = reshape_tensor(tmp, [1, [2, 3], [4, 5]])
            end
        end
    end
    return
end

"""
    function print_info(qreg::MPSQuantumRegister)

Print the dimension of virtual bond dimensions of the given MPS
"""
function print_info(qreg::MPSQuantumRegister)
    for i = 1:qreg.N
        dims = size(qreg.state[i])[2:3]
        println("$(i) has virtual dims $(dims)")
    end
end

"""
    compress!(qreg::MPSQuantumRegister)

Compress the provided MPS quantum register
"""
function compress!(qreg::MPSQuantumRegister, threshold::AbstractFloat=1e-15)
    for i in 1:qreg.N-1
        compress_bond!(qreg, i, threshold)
    end
    for i in qreg.N-1:-1:1
        compress_bond!(qreg, i, threshold)
    end
end

function compress_bond!(qreg::MPSQuantumRegister, index::Integer, threshold::AbstractFloat=1e-15)
    A = qreg.state[index]
    A_dim = size(A)
    B = qreg.state[index+1]
    B_dim = size(B)

    # Contract virtual bond connecting A and B
    @tensor C[idx1, up, idx2, down] := A[idx1, up, c] * B[idx2, c, down]

    # Reshape to a matrix and decompose
    dims = size(C)
    C = reshape(C, (dims[1]*dims[2], dims[3]*dims[4]))
    F = svd(C)

    # Take cutoff and reshape matrices back
    s = F.S
    chi = sum(s .> threshold)
    s = s[1:chi]
    qreg.s_values[index] = s./sqrt(sum(s.^2))
    A = F.U[:,1:chi] * Diagonal(sqrt.(s))
    qreg.state[index] = reshape(A, (2, A_dim[2], chi))
    B = Diagonal(sqrt.(s)) * F.Vt[1:chi,:]
    qreg.state[index+1] = permutedims(reshape(B, (chi, 2, B_dim[3])), (2, 1, 3))
end

function entropy(s::Array{<:AbstractFloat, 1})
    -sum(s.^2 .* log.(s.^2))
end

function enaglemment_entropy(qreg::MPSQuantumRegister)
    if qreg.s_values[1] == nothing
        return
    else
        return [entropy(qreg.s_values[x]) for x in 1:qreg.N-1]
    end
end

"""
    execute!(qc::QuantumCircuit, state::MPSQuantumRegister)

Specialised execute for MPS Quantum Register
"""
function execute!(qc::Gates.QuantumCircuit, state::MPSQuantumRegister,
                  compress_freq::Integer=0)
    for (i, op) in enumerate(qc.ops)
        if op.n == 1
            apply_1qubit!(state, op.gate, op.qubits[1])
        elseif op.n == 2
            apply_2qubit!(state, op.gate, op.qubits[1], op.qubits[2])
        end
        if compress_freq > 0 && i % compress_freq == 0
            compress!(state)
        end
    end
end
