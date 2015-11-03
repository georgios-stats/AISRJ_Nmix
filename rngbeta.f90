subroutine rngbeta(x,alpha,beta)

	implicit none

	double precision, intent(in)	:: alpha,beta
	double precision, intent(out)	:: x

	double precision	:: y

	call rnggamma(x,alpha)
	call rnggamma(y,beta)

	x = x/(x+y)

end subroutine rngbeta

