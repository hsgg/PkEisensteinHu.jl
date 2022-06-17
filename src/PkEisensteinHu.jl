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


module PkEisensteinHu


export compute_pk


include("growth.jl")
include("tf_fit.jl")
include("tf_fit_nu.jl")
include("eisensteinhu.jl")


function compute_pk(k_ov_h::AbstractArray, z=0.0, om0=0.3, ode0=0.7, ob0=0.05, h0=0.7, w=-1.0,
                    ns=0.9672, run=0.0, deltaR2=2e-9, mnu=93.1e-3, Nnu=3.046)
    #@show z om0 ode0 ob0 h0 w ns run deltaR2 mnu Nnu
    #@assert mnu == 0  # neutrinos not yet supported
    Tcmb = 2.726

    # set parameters
    onu0 = (mnu/93.1)/h0^2 # Omega_nu of massive neutrinos today
    fnu = onu0/om0 # neutrino fraction in mass
    pcb = (5-sqrt(25-24*fnu))/4
    zeq = 2.5e4*om0*h0^2*(2.7/2.726)^4

    # tabulate g(z) by calling "setup_growth"
    #setup_growth()
    # linear growth factor, normalized such that (1+z)D(z)=1+zeq during the matter era
    #D1=g(z)*(1+zeq)/(1+z) # now output P(k,z) as a function of k
    D1 = (1+zeq)/(1+z)

    #TFset_parameters(om0*h0^2,ob0/om0,2.726)
    tfpars = TFset_parameters_nu(om0*h0^2,ob0/om0,fnu,Nnu,Tcmb)

    om0hh = om0 * h0^2
    f_baryon = ob0 / om0
    D = D1 / (1 + zeq)

    pk = similar(k_ov_h)
    for i=1:length(pk)
        k = k_ov_h[i]
        kh0 = k * h0
        trans_nowiggle = eisensteinhu(kh0,om0hh,f_baryon) # no-wiggle P(k) with massless neutrinos
        trans_nu=TF_master(tfpars,kh0,om0hh,f_baryon,fnu,Nnu) # no-wiggle P(k) with massive neutrinos
        trans, tf_baryon, tf_cdm = TFtransfer_function(tfpars,kh0,om0hh,f_baryon) # P(k) with BAO and massless neutrinos
        trans=(trans-trans_nowiggle)+trans_nu # add BAO to no-wiggle P(k) with massive neutrinos
        # modify growth by neutrinos
        #q = k_ov_h/(om0*h0)*(2.726/2.7)^2
        #if fnu != 0   #modified by Aniket
        #    yfs = 17.2*fnu*(1+0.488/fnu^(7/6))*(Nnu*q/fnu)^2   #modified by Aniket
        #else          #modified by Aniket
        #    yfs = 1e10   #modified by Aniket
        #end         #modified by Aniket
        #D=((1-fnu)^(0.7/pcb)+(D1/(1+yfs))^0.7)^(pcb/0.7)*D1^(1-pcb)
        #D=D/(1+zeq)
        # Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.002/Mpc
        # Remember that k_ov_h is in units of h/Mpc whereas "k" in Eq.(74) is in units 
        # of 1/Mpc.
        pk[i] = (deltaR2*(2*k^2*2998^2/5/om0)^2*D^2
              *trans^2*(kh0/0.002)^(ns-1+0.5*run*log(kh0/0.002))
              *2*pi^2/k^3)
        #@show k, pk[i]
    end
    return pk
end



end


# vim: set sw=4 sts=4 et :
