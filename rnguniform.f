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
! along with AISRJ_Nmix.  If not, see <http://www.gnu.org/licenses/>.

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


      subroutine init_seedrng()
      call random_seed()
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on (0,1)-real-interval
c-----------------------------------------------------------------------
      subroutine rnguniform(rnd)
      double precision rnd
      call random_number(rnd)
      return
      end

      function runif()
      double precision runif
      call random_number(runif)
      return
      end

      function runifsp()
      real runifsp
      double precision r
      call random_number(r)
      runifsp=real(r)
      return
      end



