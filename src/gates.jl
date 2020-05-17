export Gates

module Gates

export GateOperator, QuantumCircuit
export x, y, z, h, Id, cx, cy, rk, swap

"Struct to represent operators acting on specific qubits"
mutable struct GateOperator
    n::Integer
    qubits::Array{Integer, 1}
    gate::Array{T, 2} where T <: Number
    GateOperator(qubits::Array{<:Integer, 1}, gate::Array{<:Number, 2}) =
        new(length(qubits), qubits, gate)
end

"Struct to represent a quantum circuit"
mutable struct QuantumCircuit
    n::Integer
    ops::Array{GateOperator, 1}
    QuantumCircuit(N) = new(N, Array{GateOperator, 1}())
end

"""
Define some simple gates
"""
x = [[0., 1.] [1., 0.]]
z = [[1., 0.] [0., -1.]]
y = .- z * x .* im
s = [[1., 0.] [0., 1im]]
t = [[1., 0.] [0., exp(π*1im/4)]]
h = 1/sqrt(2)*[[1., 1.] [1., -1.]]
Id = [[1., 0.] [0., 1.]]
cx = [[1., 0., 0., 0.] [0., 1., 0., 0.] [0., 0., 0., 1.] [0., 0., 1., 0.]]
cy = [[1., 0., 0., 0.] [0., 1., 0., 0.] [0., 0., 0., -1im] [0., 0., 1im, 0.]]
rk = k -> [[1., 0., 0., 0.] [0., 1., 0., 0.] [0., 0., 1., 0.] [0., 0., 0., exp(2π*1im/2^k)]]
swap = [[1. ,0., 0., 0.] [0., 0., 1., 0.] [0., 1., 0., 0] [0., 0., 0., 1.]]

# IBM physical gates
u3(θ, ϕ, λ) = [[cos(θ/2), exp(1im*λ)*sin(θ/2)] [-exp(1im*ϕ)*sin(θ/2), exp(1im*(λ + ϕ))*cos(θ/2)]]
u2(ϕ, λ) = u3(0, ϕ, λ)
u1(λ) = u3(0, 0, λ)
end # module
