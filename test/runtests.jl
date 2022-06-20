using PkEisensteinHu
using Test
using DelimitedFiles

@testset "PkEisensteinHu.jl" begin
    # compile and run test, no check for correctness, yet.
    kh = 10.0 .^ range(-5, 3, length=1000)
    z=0.0
    om0=0.306819896457
    ode0=0.693155375969
    ob0=0.048237659
    h0=0.6778
    w=-1.0
    ns=0.9672
    run=0.0
    deltaR2=1.4e-9
    mnu=93.1e-3
    Nnu=3.046
    @show om0 + ode0
    pk, pknw = compute_pk(kh; z, om0, ode0, ob0, h0, w, ns, run, deltaR2, mnu, Nnu)

    data, header = readdlm("pk.tsv"; header=true)
    kh0 = data[:,1]
    pk0 = data[:,2]
    pknw0 = data[:,3]

    @test kh0 ≈ kh  rtol=1e-15
    @test pk0 ≈ pk  rtol=1e-15
    @test pknw0 ≈ pknw  rtol=1e-15
end
