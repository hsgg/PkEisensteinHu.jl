#!/usr/bin/env julia
#   The following routines implement all of the fitting formulae in 
#   Eisenstein \& Hu (1997) 
#
#   Program TF_fit: sample driver
#   Subroutine TFset_parameters(): sets all the scalar parameters
#   Subroutine TFtransfer_function(): calculates various transfer functions
#   Functions TF_zerobaryon, TF_nowiggles, sound_horizon_fit, kpeak,
#       alpha_gamma: implement various scaling approximations of \S 4.2
#
#
# ------------------------ DRIVER ROUTINE --------------------------- 
# The following is a driver routines you might use. 
# Basically, the driver routine asks for Omega_0, the baryon fraction,
# the hubble constant, and the CMB temperature, calls TFset_parameters() to
# set all the parameters of the fit.  A loop over wavenumbers kmin to an
# inputed kmax sampled at numk per decade calls TFtransfer_function. 
#
# IMPORTANT: TFtransfer_function asks for wavenumbers in Mpc^{-1} so
#	     multiply by hubble to convert from h Mpc^{-1}
#
# The latter returns values of the various transfer functions at the given 
# wavenumber which are output to the file "trans.dat" 
#
# Also included is an example of how to call the functions
#
#	TF_nowiggles:	    shape approximation
#	TF_zerobaryon:      zero baryon TF 
#	k_peak:		    approximate first peak location
#	sound_horizon_fit:  approximate sound horizon
#	alpha_gamma:  	    small scale modification to Gamma
#
#
#
# INPUT:  omega0 -- the matter density (baryons+CDM) in units of critical 
#	  f_baryon -- the ratio of baryon density to matter density 
#	  hubble -- the Hubble constant, in units of 100 km/s/Mpc
#	  Tcmb -- the CMB temperature in Kelvin.  2.728(4) is COBE and is the default
#	          reached by setting Tcmb=0.
#	  kmax -- maximum k in  h Mpc^{-1}
#	  numk -- number of k per decade 
#
# OUTPUT: The file "trans.dat" with columns
#
#	 (1)  k -- wavenumber in h Mpc^{-1}
#	 (2)  tf_full -- The full fitting formula, eq. (16), for the matter
#			transfer function. 
#	 (3)  tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
#	 (4)  tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
#	 (5)  tf_nowiggles -- An approximate form, eqs. (30)-(31), that fits
#			only the non-oscillatory part of the transfer 
#			function.  Appropriate only for low baryon fractions.
#	 (6)  tf_zerobaryon -- The transfer function of the zero-baryon case,
#			eq. (29); i.e. what would have occured were the
#			baryons CDM instead. 			
#
#	and the approximate values of k_peak,sound_horizon,alpha_gamma 
#	to stdout, a more accurate form of sound_horizon lives in
#
#	GLOBALVARIABLES: Various intermediate fit parameters are stored in 
#                        this common block for easy access.


# ----------------------------- DRIVER ------------------------------- 


#$$$       program TFfit
#$$$
#$$$
#$$$        real    omega0,f_baryon,hubble,Tcmb,kmax,kmin
#$$$	real	k,tf_full,tf_baryon,tf_cdm,tf_nowiggles,tf_zerobaryon,
#$$$     *		k_peak
#$$$	integer numk
#$$$
#$$$c  cosmological parameters
#$$$
#$$$	write(6,*) 'Omega_0,f_baryon,h,T_cmb?'
#$$$	read*,     omega0,f_baryon,hubble,Tcmb
#$$$	omhh = omega0*hubble*hubble
#$$$
#$$$c  call routine to set fitting parameters
#$$$
#$$$        call TFset_parameters(omhh, f_baryon, Tcmb)
#$$$
#$$$
#$$$c  loop over k, call subroutine and functions to calc TFs
#$$$ 
#$$$	open(10,file='trans.dat')
#$$$
#$$$	write(6,*) 'k_max (h Mpc^{-1}),#pts per decade (10,50)'
#$$$	read*,kmax,numk
#$$$
#$$$	if (kmax.le.0) kmax=10.
#$$$	if (numk.le.0) numk=50
#$$$
#$$$	kmin = 0.0001
#$$$	numk = numk*log10(kmax/kmin)
#$$$
#$$$	do i=1,numk
#$$$
#$$$	 k=10.^(i*(log10(kmax/kmin)/numk))*kmin
#$$$         call TFtransfer_function(k*hubble,omhh,f_baryon,
#$$$     &     tf_full,tf_baryon,tf_cdm)
#$$$	 write(10,50) k,tf_full,tf_baryon,tf_cdm,
#$$$     &               TF_nowiggles(k*hubble,omhh,f_baryon,Tcmb),
#$$$     &               TF_zerobaryon(k*hubble/omhh*(Tcmb/2.7)^2)
#$$$
#$$$
#$$$	end do
#$$$
#$$$c  example of how to use the scaling functions
#$$$
#$$$	write(6,*) 'Some useful approximate scalings:'
#$$$
#$$$	write(6,10) k_peak(omhh,f_baryon)/hubble
#$$$10	FORMAT(1X,' First peak location (h Mpc^{-1}):  ',E13.4)
#$$$	write(6,20) sound_horizon_fit(omhh,f_baryon)*hubble
#$$$20	FORMAT(1X,' Approx. sound horizon (h^{-1} Mpc):',E13.4)
#$$$	write(6,30) alpha_gamma(omhh,f_baryon)
#$$$30	FORMAT(1X,' alpha_gamma:                       ',E13.4)
#$$$
#$$$
#$$$50	FORMAT(1X,7E13.5)
#$$$
#$$$      end

#
#
# PART I:------------------- FITTING FORMULAE ROUTINES ----------------- 
#
# There are two routines and a set of functions.  
#   TFset_parameters() sets all the scalar parameters, while 
#   TFtransfer_function() calculates various transfer functions 
#
# Global variables -- We've left many of the intermediate results as
# global variables in case you wish to access them, e.g. by declaring
# them as a common block in your main program. 
#
# Note that all internal scales are in Mpc, without any Hubble constants! 
#

function TFset_parameters(omhh,f_baryon,Tcmb)
	# Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
	# Input omhh -- The density of CDM and baryons, in units of critical dens,
	#                multiplied by the square of the Hubble constant, in units
	#                of 100 km/s/Mpc */
	#       f_baryon -- The fraction of baryons to CDM */
	#       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
	#		the default reached by inputing Tcmb=0 -- reset on output. */
	# Output nothing, but set many global variables in common block 
	#       GLOBALVARIABLES. You can access them yourself, if you want:
	#
	#	theta_cmb,	/* Tcmb in units of 2.7 K */ 
	#	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
	#	k_equality,	/* Scale of equality, in Mpc^-1 */
	#	z_drag,		/* Redshift of drag epoch */
	#	R_drag,		/* Photon-baryon ratio at drag epoch */
	#	R_equality,	/* Photon-baryon ratio at equality epoch */
	#	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
	#	k_silk,		/* Silk damping scale, in Mpc^-1 */
	#	alpha_c,	/* CDM suppression */
	#	beta_c,		/* CDM log shift */
	#	alpha_b,	/* Baryon suppression */
	#	beta_b,		/* Baryon envelope shift */



	# Are inputs reasonable?
	(f_baryon <= 0) && (f_baryon=1.e-5)
	(Tcmb <= 0) && (Tcmb=2.728)
        (omhh <= 0) && error("illegal input omhh=$omhh")

	#        if (hubble.gt.10.0) then
	#	   write(6,*) 'TFset_parameters(): WARNING, Hubble constant in 
	#     &	               100km/s/Mpc desired'
	#	end if

	# Auxiliary variables
        obhh = omhh*f_baryon
        theta_cmb = Tcmb/2.7

	# Main variables
        z_equality = 2.50e4*omhh*theta_cmb^(-4) - 1
        k_equality = 0.0746*omhh*theta_cmb^(-2) 

	z_drag = 0.313*omhh^(-0.419)*(1 + 0.607*omhh^(0.674))
	z_drag = 1 + z_drag*obhh^(0.238*omhh^(0.223))
        z_drag = 1291 * omhh^(0.251)/(1 + 0.659*omhh^(0.828)) * z_drag
 
        R_drag = 31.5*obhh*theta_cmb^(-4.)*1000/(1 + z_drag) 
        R_equality = 31.5*obhh*theta_cmb^(-4)*1000/(1 + z_equality) 

        sound_horizon = (2/3/k_equality*sqrt(6/R_equality)*
			 log(( sqrt(1+R_drag)+sqrt(R_drag+R_equality) )
			    /(1+sqrt(R_equality))))

        k_silk = 1.6*obhh^(0.52)*omhh^(0.73)*(1 + (10.4*omhh)^(-0.95))

	alpha_c = (46.9*omhh)^(0.670)*(1+(32.1*omhh)^(-0.532))
	alpha_c = alpha_c^(-f_baryon) 
	alpha_c = alpha_c*((12*omhh)^(0.424)*(1 + (45*omhh)^(-0.582)))^(-f_baryon^3)

    
	beta_c = 0.944/(1+(458*omhh)^(-0.708))
	beta_c = 1+beta_c*((1-f_baryon)^((0.395*omhh)^(-0.0266)) - 1)
	beta_c = 1/beta_c

	y = (1+z_equality)/(1+z_drag)
	alpha_b = y*(-6*sqrt(1+y)+(2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)))
        alpha_b = 2.07*k_equality*sound_horizon*(1+R_drag)^(-0.75)*alpha_b

        beta_b = 0.5+f_baryon+(3-2*f_baryon)*sqrt((17.2*omhh)^2+1)

        beta_node = 8.41*omhh^(0.435)

        return
end


function TFtransfer_function(tfpars, k, omhh, f_baryon)
	#  Calculate transfer function from the fitting parameters stored in
	#  GLOBALVARIABLES.
	#
	#  Input: 
	#	 k -- wavenumber in Mpc^{-1}  
	#        omhh -- The density of CDM and baryons, in units of critical dens,
	#                multiplied by the square of the Hubble constant, in units
	#                of 100 km/s/Mpc */
	#        f_baryon -- The fraction of baryons to CDM */
	#	
	#  Output:
	#	 tf_full -- The full fitting formula, eq. (16), for the matter
	#	            transfer function. 
	#	 tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
	#	 tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
	#

	#  Reasonable k?
	(k <= 0) && error("TFtransfer_function(): Illegal k")

	#  Auxiliary Variables
	q = k/13.41/tfpars.k_equality
	ks = k*tfpars.sound_horizon

	#  Main Variables
	tf_cdm = 1/(1+(ks/5.4)^4)
	tf_cdm = (tf_cdm*TF_pressureless(q,1,tfpars.beta_c)
		  + (1-tf_cdm)*TF_pressureless(q,tfpars.alpha_c,tfpars.beta_c))

	s_tilde = tfpars.sound_horizon/(1+(tfpars.beta_node/ks)^3)^(1/3) 
	tf_baryon = TF_pressureless(q,1,1)/(1+(ks/5.2)^2)
	tf_baryon = tf_baryon + tfpars.alpha_b/(1+(tfpars.beta_b/ks)^3)*exp(-(k/tfpars.k_silk)^(1.4))
	tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
	tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

	return tf_full, tf_baryon, tf_cdm
end

#       auxiliary function: Pressureless TF
function TF_pressureless(q,a,b)
	TF_pressureless = log(ℯ+1.8*b*q)
	TF_pressureless = TF_pressureless/(TF_pressureless + (14.2/a + 386/(1+69.9*q^1.08))*q^2)
	return TF_pressureless
end	
#
#
#
#
#
#
# PART II:------------------- Scaling Functions ROUTINES ----------------- 
#
#       omhh -- The density of CDM and baryons, in units of critical dens,
#                multiplied by the square of the Hubble constant, in units
#                of 100 km/s/Mpc */
#       f_baryon -- The fraction of baryons to CDM */
#
#
#	TF_zerobaryon:     
#	  Input:  q = k/omhh * (Tcmb/2.7)^2    (k in Mpc^{-1})
#	  Output: zero baryon TF Eq(29)
#	TF_nowiggles:      
#	  Input:  k = wavenumber in Mpc^{-1}, omhh, f_baryon, Tcmb
#	  Output: shape approximation TF  Eq(30-31)
#	  Calls: TF_zerobaryon,sound_horizon_fit,alpha_gamma
# 	sound_horizon_fit: 
#         Input:  omhh,f_baryon	
#	  Output: approximate sound horizon in Mpc	
#	kpeak:		   
#	  Input:  omhh,f_baryon
#         Output: first peak location in Mpc
#	  Calls:  sound_horizon_fit
#	alpha_gamma:	   
#	  Input: omhh,f_baryon
#	  Output: effective small scale suppression


function TF_zerobaryon(q)
	TF_zerobaryon = log(2*e+1.8*q)
	TF_zerobaryon = TF_zerobaryon/(TF_zerobaryon+(14.2 + 731.0/(1+62.5*q))*q^2)
	return
end

function TF_nowiggles(k,omhh,f_baryon,Tcmb)
	(Tcmb <= 0) && (Tcmb=2.728)
	a = alpha_gamma(omhh,f_baryon)
	q_eff = k/omhh*(Tcmb/2.7)^2
	q_eff = q_eff/(a+(1-a)/(1+(0.43*k*sound_horizon_fit(omhh,f_baryon))^4))
	TF_nowiggles = TF_zerobaryon(q_eff)
	return
end


function sound_horizon_fit(omhh,f_baryon)
	obhh = f_baryon*omhh
	sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*obhh^(0.75))
	return
end


function k_peak(omhh,f_baryon)
	obhh = f_baryon*omhh
	k_peak = 5*pi/2*(1+0.217*omhh)/sound_horizon_fit(omhh,f_baryon)
	return
end


function alpha_gamma(omhh,f_baryon)
	alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*(f_baryon)^2
	return
end 
