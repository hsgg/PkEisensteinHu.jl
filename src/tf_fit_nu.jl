#!/usr/bin/env julia
#	Subroutines extracted from Wayne Hu's "power.f" available at
#     http://background.uchicago.edu/~whu/transfer/power.f

struct TFEisensteinHu_nu
	theta_cmb::Float64
	alpha_nu::Float64
	alpha_b::Float64
	beta_b::Float64
	alpha_c::Float64
	beta_c::Float64
	beta_node::Float64
	sound_horizon::Float64
	k_equality::Float64
	k_silk::Float64
end

function TFset_parameters_nu(omhh,f_baryon,f_nu,N_nu,Tcmb)
	# Auxiliary variable
	obhh = omhh*f_baryon
	theta_cmb = Tcmb/2.7

	# Main variables
	z_equality = 2.50e4*omhh*theta_cmb^(-4) - 1
	k_equality = 0.0746*omhh*theta_cmb^(-2)

	z_drag = 0.313*omhh^(-0.419)*(1+0.607*omhh^(0.674))
	z_drag = 1 + z_drag*obhh^(0.238*omhh^(0.223))
	z_drag = 1291 * omhh^(0.251)/(1 + 0.659*omhh^(0.828)) * z_drag

	y_d = (1+z_equality)/(1+z_drag)

	R_drag = 31.5*obhh*theta_cmb^(-4)*1000e0/(1 + z_drag)
	R_equality = 31.5*obhh*theta_cmb^(-4)*1000/(1 + z_equality)

	sound_horizon = (2/3/k_equality*sqrt(6/R_equality)*
			 log(( sqrt(1+R_drag)+sqrt(R_drag+R_equality) )
			    /(1+sqrt(R_equality))))

	p_c  = -(5-sqrt(1+24*(1-f_nu-f_baryon)))/4
	p_cb = -(5-sqrt(1+24*(1-f_nu)))/4
	f_c  = 1-f_nu-f_baryon
	f_cb = 1-f_nu
	f_nub= f_nu+f_baryon

	alpha_nu = (f_c/f_cb) * (2*(p_c+p_cb)+5)/(4*p_cb+5)
	alpha_nu = alpha_nu*(1-0.553*f_nub+0.126*f_nub^3)
	alpha_nu = alpha_nu/(1-0.193*sqrt(f_nu)+0.169*f_nu)
	alpha_nu = alpha_nu*(1+y_d)^(p_c-p_cb)
	alpha_nu = alpha_nu*(1+ (p_cb-p_c)/2 * (1+1/(4*p_c+3)/(4*p_cb+7))/(1+y_d))

	y = (1+z_equality)/(1+z_drag)
	alpha_b = y*(-6*sqrt(1+y)+(2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)))
        alpha_b = 2.07*k_equality*sound_horizon*(1+R_drag)^(-0.75)*alpha_b

        beta_b = 0.5+f_baryon+(3-2*f_baryon)*sqrt((17.2*omhh)^2+1)

	alpha_c = (46.9*omhh)^(0.670)*(1+(32.1*omhh)^(-0.532))
	alpha_c = alpha_c^(-f_baryon) 
	alpha_c = alpha_c*((12*omhh)^(0.424)*(1 + (45*omhh)^(-0.582)))^(-f_baryon^3)

	beta_c=1/(1-0.949*f_nub)

        beta_node = 8.41*omhh^(0.435)

        k_silk = 1.6*obhh^(0.52)*omhh^(0.73)*(1 + (10.4*omhh)^(-0.95))

	return TFEisensteinHu_nu(theta_cmb, alpha_nu, alpha_b, beta_b, alpha_c, beta_c, beta_node, sound_horizon, k_equality, k_silk)
end 


function TF_master(tfpars,k,omhh,f_baryon,f_nu,N_nu)
	q = k*tfpars.theta_cmb^2/omhh
	gamma_eff=sqrt(tfpars.alpha_nu) + (1-sqrt(tfpars.alpha_nu))/(1+(0.43*k*tfpars.sound_horizon)^4)

	q_eff = q/gamma_eff
	TF_master= log(ℯ+1.84*tfpars.beta_c*sqrt(tfpars.alpha_nu)*q_eff)
	TF_master = TF_master/(TF_master + q_eff^2*(14.4 + 325/(1+60.5*q_eff^1.11)))

	q_nu = 3.92*q*sqrt(N_nu/f_nu)
	TF_master = TF_master*(1+(1.2*f_nu^(0.64)*N_nu^(0.3+0.6*f_nu))/(q_nu^(-1.6)+q_nu^(0.8)))
	return TF_master
end
