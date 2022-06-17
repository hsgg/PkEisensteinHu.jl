#!/usr/bin/env julia


module TestPlotPk

using PkEisensteinHu
using PyPlot
using DelimitedFiles


function plot_pk()
    plt.close("all")

    k = 10.0 .^ range(-5, 3, length=1000)
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
    #mnu=0.0
    Nnu=3.046
    @show om0 + ode0
    pk0, pk1 = compute_pk(k, z, om0, ode0, ob0, h0, w, ns, run, deltaR2, mnu, Nnu)

    open("pk.tsv", "w") do f
        write(f, "# kh    pk  pk_nowiggle\n")
        writedlm(f, [k pk0 pk1])
    end


    mask = @. (1e-4 <= k <= 1)
    k = k[mask]
    pk0 = pk0[mask]
    pk1 = pk1[mask]


    figure()
    plot(k, pk0, c="k", label="fiducial")
    plot(k, pk1, c="C0", label="nowiggle")
    for i=-10:20:10
        h = h0 + i / 100
        pk, _ = compute_pk(k, z, om0, ode0, ob0, h, w, ns, run, deltaR2, mnu, Nnu)
        plot(k, pk, c=string(h), alpha=0.5)
    end
    legend()
    xlabel(L"k")
    ylabel(L"P(k)")
    xscale("log")
    yscale("log")


    figure()
    plot(k, pk0 .- pk1, c="k", label="diff")
    #plot(k, 100*sin.(k .* 151 * 0.7))
    legend()
    xlabel(L"k")
    ylabel(L"P(k) - P_{\rm no\,wiggle}(k)")
    #xscale("log")
    #yscale("log")


    figure()
    plot(k, pk0 ./ pk1, c="k", label="ratio")
    legend()
    xlabel(L"k")
    ylabel(L"P(k) / P_{\rm no\,wiggle}(k)")
    #xscale("log")
    #yscale("log")
end


end

TestPlotPk.plot_pk()


# vim: set sw=4 et sts=4 :
