
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

@testset "Test MPS compression on full MPS register" begin
    @test begin
        qreg = MPSQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.h, 1)
        apply_2qubit!(qreg, Gates.cx, 1, 2)
        apply_2qubit!(qreg, Gates.cx, 2, 3)
        qreg_state = convert(FullStateQuantumRegister{ComplexF64}, qreg)
        compress!(qreg)
        convert(FullStateQuantumRegister{ComplexF64}, qreg) ≈ qreg_state
    end

    @testset "On larger system with random prep circuit" begin
        qreg_mps = MPSQuantumRegister{ComplexF64}(10, "0"^10)
        qreg_full = FullStateQuantumRegister{ComplexF64}(10, "0"^10)
        qc = simple_state_prep_circuit(10, 8)
        execute!(qc, qreg_mps)
        execute!(qc, qreg_full)

        @test begin
            # check that before compression the full state equals mps state
            qreg_full ≈ convert(FullStateQuantumRegister{ComplexF64}, qreg_mps)
        end

        @test begin
            # check states still equal after compression
            compress!(qreg_mps)
            qreg_full ≈ convert(FullStateQuantumRegister{ComplexF64}, qreg_mps)
        end
    end
end

@testset "Test execute function with compression" begin
    qreg1 = MPSQuantumRegister{ComplexF64}(3, "000")
    qreg2 = MPSQuantumRegister{ComplexF64}(3, "000")
    qc = ghz_circuit(3)
    # add some additional unecessary gates that should cancel out
    push!(qc.ops, Gates.GateOperator([1, 2], Gates.cx))
    push!(qc.ops, Gates.GateOperator([1, 2], Gates.cx))
    push!(qc.ops, Gates.GateOperator([2, 3], Gates.cx))
    push!(qc.ops, Gates.GateOperator([2, 3], Gates.cx))
    push!(qc.ops, Gates.GateOperator([1, 3], Gates.cx))
    push!(qc.ops, Gates.GateOperator([1, 3], Gates.cx))    
    execute!(qc, qreg1) # with no compression
    execute!(qc, qreg2, 1) # with compression after each gate application

    @test begin
        # test two registers match
        convert(FullStateQuantumRegister{ComplexF64}, qreg1) ≈
        convert(FullStateQuantumRegister{ComplexF64}, qreg2)
    end

    @test begin
        # test that bond dimenions are smaller
        all(size(qreg1.state[2])[2:end] .> size(qreg2.state[2])[2:end])
    end
end
