#!/usr/bin/env julia


module TestPlotPk

using PkEisensteinHu
using PyPlot


function plot_pk()
    plt.close("all")

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
    pkeh0 = compute_pk(k, z, om0, ode0, ob0, h0, w, ns, run, deltaR2, mnu, Nnu)

    figure()
    plot(k, pkeh0, c="k", label="fiducial")
    for i=-10:10
        h = h0 + i / 100
        pkeh = compute_pk(k, z, om0, ode0, ob0, h, w, ns, run, deltaR2, mnu, Nnu)
        plot(k, pkeh, c=string(h), alpha=0.5)
    end
    legend()
    xlabel(L"k")
    ylabel(L"P(k)")
    xscale("log")
    yscale("log")
end


end

TestPlotPk.plot_pk()


# vim: set sw=4 et sts=4 :
