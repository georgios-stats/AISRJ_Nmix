! CONTACT DETAILS :

! Georgios Karagiannis © 2012
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email (current): Georgios.Karagiannis@pnnl.gov

! Christophe Andrieu
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email: C.Andrieu@bristol.ac.uk

! ----------------------------------------------------------------------

subroutine rngbeta(x,alpha,beta)

	implicit none

	double precision, intent(in)	:: alpha,beta
	double precision, intent(out)	:: x

	double precision	:: y

	call rnggamma(x,alpha)
	call rnggamma(y,beta)

	x = x/(x+y)

end subroutine rngbeta

