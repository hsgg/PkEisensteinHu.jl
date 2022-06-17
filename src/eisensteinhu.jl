function eisensteinhu(ak,omegamh2,fb)
	# REF: Eisenstein and Hu, ApJ, 496, 605 (1998), Eq.(29)-(31)
	alpha= 1 - 0.328*log(431*omegamh2)*fb + 0.38*log(22.3*omegamh2)*fb^2
	sound= 44.5*log(9.83/omegamh2)/sqrt(1+10*(fb*omegamh2)^(3/4))
	shape= omegamh2*(alpha+(1-alpha)/(1+(0.43*ak*sound)^4)) 
	aq= ak*(2.725/2.7)^2/shape
	T= log(2*ℯ+1.8*aq)/(log(2*ℯ+1.8*aq)+(14.2+731/(1+62.5*aq))*aq*aq)
	return T
end
