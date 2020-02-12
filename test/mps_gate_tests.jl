
@testset "Test MPS register initialisation" begin
    @test begin
        qreg = MPSQuantumRegister{Int}(3, "000")
        to_str(qreg) == "(1|000>)"
    end
    @test begin
        qreg = MPSQuantumRegister{Int}(3, "100")
        to_str(qreg) == "(1|100>)"
    end
    @test begin
        qreg = MPSQuantumRegister{Int}(3, "001")
        to_str(qreg) == "(1|001>)"
    end
end

@testset "Test convert MPS register to FullStateRegister" begin
    @test begin
        qreg = MPSQuantumRegister{Int}(3, "100")
        qreg_ref = FullStateQuantumRegister{Int}(3, "100")
        convert(FullStateQuantumRegister{Int}, qreg) == qreg_ref
    end
    @test begin
        qreg = MPSQuantumRegister{Int}(3, "100")
        qreg_ref = FullStateQuantumRegister{ComplexF64}(3, "100")
        convert(FullStateQuantumRegister{ComplexF64}, qreg) == qreg_ref
    end
end

@testset "Test applying 1 qubit gates on mps quantum register" begin
    @test begin
        qreg = MPSQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.x, 3)
        to_str(qreg) == to_str(FullStateQuantumRegister{ComplexF64}(3, "001"))
    end
end

@testset "Test conversion from two qubit gate to mpo and back" begin
    @test begin
        mpo = two_qubit_gate_to_mpo(Gates.cx)
        gate = mpo_to_two_qubit_gate(mpo)
        gate ≈ Gates.cx
    end
    @test begin
        mpo = two_qubit_gate_to_mpo(Gates.rk(4))
        gate = mpo_to_two_qubit_gate(mpo)
        gate ≈ Gates.rk(4)
    end
    @test begin
        gate = rand(ComplexF64, 4, 4)
        gate = gate + gate'
        mpo = two_qubit_gate_to_mpo(gate)
        gate_from_mpo = mpo_to_two_qubit_gate(mpo)
        gate ≈ gate_from_mpo
    end
end

@testset "Test applying 2 qubit gates on mps quantum register" begin
    @test begin
        qreg = MPSQuantumRegister{ComplexF64}(3, "100")
        apply_2qubit!(qreg, Gates.cx, 1, 2)
        convert(FullStateQuantumRegister{ComplexF64}, qreg) ≈ FullStateQuantumRegister{ComplexF64}(3, "110")
    end
    @test begin
        qreg = MPSQuantumRegister{ComplexF64}(3, "100")
        apply_2qubit!(qreg, Gates.cx, 2, 1)
        convert(FullStateQuantumRegister{ComplexF64}, qreg) ≈ FullStateQuantumRegister{ComplexF64}(3, "100")
    end
end
