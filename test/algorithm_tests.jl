using FFTW

@testset "Algorithm tests" begin
    @test begin
        n = 2
        results = zeros(Bool, 2^n)
        for i in 1:2^n
            in_state = zeros(Float64, 2^n)
            in_state[i] = 1.0
            qreg = FullStateQuantumRegister{ComplexF64}(2, in_state)
            qc = qft_circuit(n)
            execute!(qc, qreg)
            ref_state = ifft(in_state)
            ref_state /= sqrt(sum(conj(ref_state) .* ref_state))
            qreg_ref = FullStateQuantumRegister{ComplexF64}(n, ref_state)
            results[i] = qreg ≈ qreg_ref
        end
        all(results)
    end
    @test begin
        n = 4
        in_state = rand(2^n)
        in_state /= sqrt(sum(in_state .* in_state))
        qreg = FullStateQuantumRegister{ComplexF64}(n, in_state)
        qc = qft_circuit(n)
        execute!(qc, qreg)
        ref_state = ifft(in_state)
        ref_state /= sqrt(sum(conj(ref_state) .* ref_state))
        qreg_ref = FullStateQuantumRegister{ComplexF64}(n, ref_state)
        qreg ≈ qreg_ref
    end
    @test begin
        n = 10
        in_state = rand(2^n)
        in_state /= sqrt(sum(in_state .* in_state))
        qreg = FullStateQuantumRegister{ComplexF64}(n, in_state)
        qc = qft_circuit(n)
        execute!(qc, qreg)
        ref_state = ifft(in_state)
        ref_state /= sqrt(sum(conj(ref_state) .* ref_state))
        qreg_ref = FullStateQuantumRegister{ComplexF64}(n, ref_state)
        qreg ≈ qreg_ref
    end
end
