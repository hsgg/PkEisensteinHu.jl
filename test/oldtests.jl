#!/usr/bin/env julia
# A sample program for computing the linear power spectrum
# of density fluctuations using Eisenstein&Hu's transfer function
# that includes the effects of massive neutrinos and
# the baryonic oscillation, multiplied by the growth function
# that includes the effect of massive neutrinos. 
# The modification due to massive neutrinos is given by subroutines 
# extracted from Wayne Hu's "power.f" (given in tf_fit_nu.f), and 
# the growth  function fit derived in Hu & Eisenstein (1998).
# The baryon acoustic oscillation is included by computing the difference
# with respect to the no-wiggle power spectrum for massless neutrinos.
# Ref: Eisenstein & Hu, ApJ, 496, 605 (1998) for transfer functions with massless nu
#      Eisenstein & Hu, ApJ, 511, 5 (1999) for no-wiggle transfer function with massive nu
# - k is in units of h Mpc^-1
# - P(k) is in units of h^-3 Mpc^3
# November 23, 2015: E.Komatsu
# March 4, 2016: modified by Aniket Agrawal to output the scale-dependent growth factors
# March 13, 2018: ported from FORTRAN to Julia 0.6 by Henry Gebhardt, for fnu=0

include("growth.jl")
include("../../myjl/include.jl")

module ComputePk

include("tf_fit.jl")
include("tf_fit_nu.jl")
include("eisensteinhu.jl")

using PyPlot
using AmoebaOptim
using DelimitedFiles


function test_compute_pk()
    # CAMB
    tmp = readdlm("planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
    k = tmp[:,1]
    pkcamb = tmp[:,2]

    # Eisenstein+Hu
    #d, h = readdlm("cosmo2/cosmo.tsv", header=true)
    #om0 = d[1,end-1]
    #ode0 = d[1,end-3]
    #ob0 = 0.05
    #H0 = 100e5 / 3e10
    #h0 = 0.6778
    #@show om0 ode0 ob0 h0
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
    pkeh = compute_pk(k, z, om0, ode0, ob0, h0, w, ns, run, deltaR2, mnu, Nnu)
    @show pkeh

    writedlm("cosmo2/pkeh.tsv", [k pkeh])

    figure()
    plot(k, pkcamb, label="CAMB")
    plot(k, pkeh, label="Eisenstein+Hu")
    legend()
    xlabel(L"k")
    ylabel(L"P(k)")
    xscale("log")
    yscale("log")
end

function fit_compute_pk()
    # CAMB
    tmp = readdlm("planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
    k = tmp[:,1]
    pkcamb = tmp[:,2]
    lnpkcamb = log.(pkcamb)

    # Eisenstein+Hu
    d, h = readdlm("cosmo2/cosmo.tsv", header=true)
    om0 = d[1,end-1]
    ode0 = d[1,end-3]
    ob0 = 0.05
    H0 = 100e5 / 3e10
    h0 = 0.6778
    @show om0 ode0 ob0 h0

    pfull = [0.0, om0, ode0, ob0, h0, -1.0, 0.9672, 0.0, 1.4e-9, 93.1e-3, 3.046]
    varpars = [2, 4, 5, 7, 9]
    nlnLfn(p) = begin
        pfull[varpars] .= p
        pfull[3] = 1 - pfull[2]
        (0.1 < pfull[2] < 0.4) || return 1e10
        (0.5 < pfull[3] < 0.9) || return 1e10
        (0.001 < pfull[4] < 0.1) || return 1e10
        (0.6 < pfull[5] < 0.8) || return 1e10
        (0.9 < pfull[7] < 1.1) || return 1e10
        (0.1e-9 < pfull[9] < 10e-9) || return 1e10
        pkeh = compute_pk(k, pfull...)
        fom = @. (log(pkeh) - lnpkcamb)^2
        r = sum(fom)
        @show pfull, r
        return r
    end

    dpfull = 0.01 * pfull
    dpfull[8] = 0.01
    dp = dpfull[varpars]
    pinit = collect(pfull[varpars])
    @time pmin = amoeba(pinit, dp, nlnLfn, nmax=100_000, ftol=100eps(1.0))

    pfullmin = copy(pfull)
    pfullmin[varpars] .= pmin
    pkeh = compute_pk(k, pfullmin...)
    @show pfullmin pkeh

    figure()
    plot(k, pkcamb, label="CAMB")
    plot(k, pkeh, label="Eisenstein+Hu")
    legend()
    xlabel(L"k")
    ylabel(L"P(k)")
    xscale("log")
    yscale("log")
end


end  # module ComputePk

# vim: set sw=4 et sts=4 :
