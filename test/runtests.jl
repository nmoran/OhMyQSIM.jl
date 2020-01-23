using OhMyQSIM
using Test
tests = ["gate_tests", "algorithm_tests"]
for t in tests
    include("$(t).jl")
end
