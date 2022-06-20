using PkEisensteinHu
using Test

@testset "PkEisensteinHu.jl" begin
    # compile and run test, no check for correctness, yet.
    k = 10.0 .^ range(-5, 1, length=1000)
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
    pkeh0 = compute_pk(k; z, om0, ode0, ob0, h0, w, ns, run, deltaR2, mnu, Nnu)
end
