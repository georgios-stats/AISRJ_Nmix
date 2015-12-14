

! Copyrigtht 2012 Georgios Karagiannis
!
! This file is part of AISRJ_Nmix.
!
! AISRJ_Nmix is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation version 2 of the License.
!
! AISRJ_Nmix is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

! ----------------------------------------------------------------------

! REFERENCES:
! -----------
!
! Karagiannis, G., & Andrieu, C. (2013).
! Annealed importance sampling reversible jump MCMC algorithms.
! Journal of Computational and Graphical Statistics, 22(3), 623-648.
!
! Georgios Karagiannis
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email : Georgios.Karagiannis@pnnl.gov
! Email (current): georgios-stats@gmail.com
!
! Christophe Andrieu
! School of Mathematics, University of Bristol
! University Walk, Bristol, BS8 1TW, UK
! Email: C.Andrieu@bristol.ac.uk

subroutine rngbeta(x,alpha,beta)

	implicit none

	double precision, intent(in)	:: alpha,beta
	double precision, intent(out)	:: x

	double precision	:: y

	call rnggamma(x,alpha)
	call rnggamma(y,beta)

	x = x/(x+y)

end subroutine rngbeta

