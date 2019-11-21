@testset begin
    @test all(FullStateQuantumRegister{ComplexF64}(3, "000").state .=== convert(Array{ComplexF64, 1}, [1, 0, 0, 0, 0, 0, 0, 0]))
    @test all(FullStateQuantumRegister{ComplexF64}(3, "111").state .=== convert(Array{ComplexF64, 1}, [0, 0, 0, 0, 0, 0, 0, 1]))
    @test all(FullStateQuantumRegister{ComplexF16}(3, "101").state .=== convert(Array{ComplexF16, 1}, [0, 0, 0, 0, 0, 1, 0, 0]))
    @test begin
        qreg = FullStateQuantumRegister{Int}(3, "000")
        to_str(qreg) == "(1|000>)"
    end
    @test begin
        qreg = FullStateQuantumRegister{Int}(3, "000")
        qreg.state[end] = 1
        to_str(qreg) == "(1|000>) + (1|111>)"
    end
    @test begin
        a = FullStateQuantumRegister{ComplexF64}(3, "000")
        b = FullStateQuantumRegister{ComplexF64}(3, "000")
        b.state *= 1im
        a ≈ b
    end

    @test begin
        a = FullStateQuantumRegister{ComplexF64}(3, "000")
        b = FullStateQuantumRegister{ComplexF64}(3, "000")
        b.state *= 1im
        (a == b) == false
    end

end


@testset begin
    @test begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit(qreg, Gates.x, 1) == FullStateQuantumRegister{ComplexF64}(3, "100").state
    end
    @test begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.x, 1)
        apply_1qubit(qreg, Gates.z, 1) == -FullStateQuantumRegister{ComplexF64}(3, "100").state
    end
    @test begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.y, 1)
        apply_1qubit!(qreg, Gates.x, 2)
        apply_1qubit(qreg, Gates.y, 2) == FullStateQuantumRegister{ComplexF64}(3, "100").state
    end
    @testset begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        random_gate = rand(ComplexF64, (2, 2))
        for i = 1:3
            @test apply_1qubit_full(qreg, random_gate, i) == apply_1qubit(qreg, random_gate, i)
        end
    end
end

# Test by decomposing random matrices and recomposing
@test begin
    rmat = rand(ComplexF64, (4, 4))
    ops = decompose_2_qubit_gate(rmat)
    rrmat = nothing
    for op in ops
        mat = kron(op[2], op[1])
        if rrmat == nothing
            rrmat = mat
        else
            rrmat += mat
        end
    end
    rrmat ≈ rmat
end

@test kron(Gates.Id, Gates.x) == swap_2qubits(kron(Gates.x, Gates.Id))
@test begin
    A = rand(ComplexF64, (2, 2))
    B = rand(ComplexF64, (2, 2))
    kron(A, B) == swap_2qubits(kron(B, A))
end

@testset begin
    @test begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.x, 1)
        apply_2qubit!(qreg, Gates.cx, 2, 3)
        apply_2qubit!(qreg, Gates.cx, 1, 2)
        qreg ≈ FullStateQuantumRegister{ComplexF64}(3, "110")
    end

    @test begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.h, 1)
        apply_2qubit_full(qreg, Gates.cy, 1, 2) ≈ apply_2qubit(qreg, Gates.cy, 1, 2)
    end

end

@testset begin
    @testset begin
        qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
        apply_1qubit!(qreg, Gates.h, 1)
        apply_2qubit!(qreg, Gates.cx, 1, 2)
        apply_2qubit!(qreg, Gates.cx, 1, 3)
        @test measure_probs(qreg) ≈ [0.5, 0., 0., 0., 0., 0., 0., 0.5]
        counts = get_counts(measure(qreg, shots=100))
        @test counts["000"] + counts["111"] == 100
    end

    @testset begin
        cprobs = [0.5, 0.6, 0.75, 1.0]
        @test get_conf(cprobs, 0, 2) == "00"
        @test get_conf(cprobs, 0.25, 2) == "00"
        @test get_conf(cprobs, 0.5, 2) == "00"
        @test get_conf(cprobs, 0.65, 2) == "10"
        @test get_conf(cprobs, 0.99, 2) == "11"
    end


end
