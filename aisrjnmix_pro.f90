

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


! ----------------------------------------------------------------------

program aisrjnmix_pro

! ADD THE MODULES ======================================================

      use fixedmove_mod, only             : FixedMoveSweep
      use splitmergemove_mod, only        : logImportanceWeightSplit, &
                                          AisSplitSweep

! DECLARE THE PARAMETERS ===============================================

      implicit none

      ! FOR THE NUMBER GENERATOR

      integer                 :: seed

      ! PARAMETERS OF THE ALGORITHM (GENERAL)

      integer                 :: nsweeps,nburnin
      integer                 :: enmax,kmaxmax
      integer                 :: kmax,kmin
      parameter               ( enmax = 500 )
      parameter               ( kmaxmax = 50 )

      ! PARAMETERS OF THE ALGORITHM (AIS)

      integer                 :: Tau

      double precision        :: scaleW
      double precision        :: scaleM
      double precision        :: scaleS

      double precision        :: scaleWais
      double precision        :: scaleMais
      double precision        :: scaleSais

      ! DATASET

      character(len=80)       :: data_path
      integer                 :: en
      double precision        :: y(enmax)

      ! RANDOM PARAMETERS

      integer                 :: k

      double precision        :: w(kmaxmax)
      double precision        :: m(kmaxmax)
      double precision        :: s(kmaxmax)
      double precision        :: beta

      ! RANDOM PARAMETERS AIS
      integer                 :: kais
      double precision        :: wais(kmaxmax)
      double precision        :: mais(kmaxmax)
      double precision        :: sais(kmaxmax)
      double precision        :: betaais

      ! AUXILIARY VARIABLES
      double precision        :: wj0,mj0,sj0
      double precision        :: wj1,mj1,sj1
      double precision        :: wj2,mj2,sj2
      double precision        :: waux,maux,saux
      logical                 :: qadjacent

      ! FIXED PARAMETERS
      double precision        :: delta
      double precision        :: xi
      double precision        :: kappa
      double precision        :: alpha
      double precision        :: ge
      double precision        :: ha

      ! RJ PROPOSALS
      double precision        :: qsplit(kmaxmax)
      double precision        :: qmerge(kmaxmax)

      ! AIS SWEEP PARAMETERS
      integer                 :: t
      double precision        :: gt,Dgt
      double precision        :: logAIW,logIWt

      ! ESTIMATIONS

      integer                 :: nSplit
      integer                 :: nMerge
      integer                 :: nSplitNadj

      integer                 :: sumK
      integer                 :: sumKsq

      double precision        :: ExpAccPrSplit
      double precision        :: ExpAccPrMerge

      double precision        :: ExpAccPrW
      double precision        :: ExpAccPrM
      double precision        :: ExpAccPrS

      double precision        :: ExpAccPrWais
      double precision        :: ExpAccPrMais
      double precision        :: ExpAccPrSais

      ! LABELS

      integer                 :: k1
      integer                 :: k2

      ! SEEDS

      !OTHERS
      integer                 :: i,j
      integer                 :: iter
      integer                 :: idisp
      double precision        :: ymax,ymin
      double precision        :: u
      integer                 :: j0

      double precision        :: AccPrRJ

      double precision        :: AccPrW
      double precision        :: AccPrM
      double precision        :: AccPrS

      double precision        :: AccPrWais
      double precision        :: AccPrMais
      double precision        :: AccPrSais

      ! ARGUMENTS
      integer                 :: ia,iargc
      character(len=80)       :: word
      integer                 :: istdTau
      character(len=14)       :: specTau

! PARAMETERS ===========================================================

      ! DEFAULT ONES

      nsweeps = 2*10**5
      nburnin = max(10**4,nsweeps/10)

      Tau = 1

      data_path = "./dataset/enzymedata"

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:5).eq.'-Tau=') then
                  read(word,'(5x,i14)') Tau
            else if(word(1:9).eq.'-Nsweeps=') then
                  read(word,'(9x,i14)') nsweeps
                  nburnin = max(10**4,nsweeps/10)
            else if(word(1:9).eq.'-Nburnin=') then
                  read(word,'(9x,i14)') nburnin
            else if(word(1:1).ne.'-') then
                  data_path = word
            end if
      end do

! OPEN THE FILES =======================================================

      istdTau = 14-int(log10(0.5d0+Tau))
      write(specTau,'(i14)') Tau

      ! LOG FILES
      open(1, file = './results/log'//'.T'//specTau(istdTau:14))
      ! DATASETS
      open(2, file = data_path,status='old')
      ! MODEL LABEL
      open(3, file = './results/k'//'.T'//specTau(istdTau:14))
      ! MODEL PARAMETERS
      open(4, file = './results/w'//'.T'//specTau(istdTau:14))
      open(5, file = './results/m'//'.T'//specTau(istdTau:14))
      open(6, file = './results/s'//'.T'//specTau(istdTau:14))
      open(7, file = './results/b'//'.T'//specTau(istdTau:14))
      ! LOCAL Split/Merge MOVE
      open(8, file = './results/move'//'.T'//specTau(istdTau:14))
      open(9, file = './results/logAIW'//'.T'//specTau(istdTau:14))
      open(10, file = './results/AccPrRJ'//'.T'//specTau(istdTau:14))

      ! REPORT LOG -----------------------------------------------------

      do idisp = 0,1
      write(idisp,*) '   '
      write(idisp,*) 'ALGORITHM PARAMETERS ============================'
      write(idisp,*) '   '
      write(idisp,*) 'BURN IN AREA                    : ', nburnin
      write(idisp,*) '   '
      write(idisp,*) 'NUMBER OF ITERATIONS            : ', nsweeps
      write(idisp,*) '   '
      write(idisp,*) 'Tau                             : ', Tau
      write(idisp,*) '   '
      end do

! START THE RANDOM NUMBER GENERATOR ====================================

      call system_clock(count=seed)
      seed = seed/2
      call init_genrand(seed)
      do i = 1,5
            call rnguniform(u)
      end do

      ! REPORT LOG -----------------------------------------------------

      do idisp = 0,1
      write(idisp,*) ' '
      write(idisp,*) 'RNG SEED                        : ', seed
      write(idisp,*) ' '
      end do

! READ THE DATASET =====================================================

      read(2,*) en
      do i = 1, en
            read(2,'(f18.14)') y(i)
      end do

      ! PRINT DATAPATH & DATASET
      do idisp = 0,1
      write(idisp,*) ' '
      write(idisp,*) 'DATAPATH                        : ', data_path
      write(idisp,*) ' '
      write(idisp,'(500f16.9)') y(1:en)
      write(idisp,*) ' '
      end do

! COMPUTE THE FIXED PARAMETERS =========================================

      ymax = y(1)
      ymin = y(1)
      do i = 2,en
            if (y(i).gt.ymax) ymax = y(i)
            if (y(i).lt.ymin) ymin = y(i)
      end do
      kmax = 30
      kmin = 1
      xi = 0.5d0*(ymax+ymin)
      kappa = 1.0d0/(ymax-ymin)**2
      alpha = 2.0d0
      ge = 0.2d0
      ha = 10.d0/(ymax-ymin)**2
      delta = 1.0d0

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:6).eq.'-kmin=') then
                  read(word,'(6x,i14)') kmin
            else if(word(1:6).eq.'-kmax=') then
                  read(word,'(6x,i14)') kmax
            else if(word(1:4).eq.'-xi=') then
                  read(word,'(4x,f16.10)') xi
            else if(word(1:7).eq.'-kappa=') then
                  read(word,'(7x,f16.10)') kappa
            else if(word(1:7).eq.'-alpha=') then
                  read(word,'(7x,f16.10)') alpha
            else if(word(1:4).eq.'-ge=') then
                  read(word,'(4x,f16.10)') ge
            else if(word(1:4).eq.'-ha=') then
                  read(word,'(4x,f16.10)') ha
            else if(word(1:7).eq.'-delta=') then
                  read(word,'(7x,f16.10)') delta
            end if
      end do

      ! REPORT LOG -----------------------------------------------------

      do idisp = 0,1
      write(idisp,*) '   '
      write(idisp,*) 'FIXED PARAMETERS (MODEL) ========================'
      write(idisp,*) 'Kmax                            : ', kmax
      write(idisp,*) 'Kmin                            : ', kmin
      write(idisp,*) 'DELTA                           : ', delta
      write(idisp,*) 'XI                              : ', xi
      write(idisp,*) 'KAPPA                           : ', kappa
      write(idisp,*) 'ALPHA                           : ', alpha
      write(idisp,*) 'GE                              : ', ge
      write(idisp,*) 'HA                              : ', ha
      write(idisp,*) '   '
      end do

! COMPUTE THE RJ PROPOSALS =============================================

      do j = kmin+1,kmax-1
            qsplit(j) = 0.5d0
            qmerge(j) = 0.5d0
      end do
      qsplit(kmin) = 1.0d0
      qsplit(kmax) = 0.0d0
      qmerge(kmin) = 0.0d0
      qmerge(kmax) = 1.0d0

      ! REPORT LOG -----------------------------------------------------

      do idisp = 0,1
      write(idisp,*) ' '
      write(idisp,*) 'RJ PROPOSALS ===================================='
      write(idisp,*) '       ',' ','      Split/Merge    '
      write(idisp,*) 'k      ',' ','P(k->k+1) ',' ','P(k->k-1) '
      do k = kmin,kmax
      write(idisp,'(i2,4X,6(f10.5,X))') k,qsplit(k),qmerge(k)
      end do
      write(idisp,*) ' '
      end do

! SCALE THE FIXED DIMENSION UPDATES (GIBBS) ============================

      ! DEFAULTS
      scaleW = 0.27d0
      scaleM = 0.018d0
      scaleS = 0.29d0

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:8).eq.'-scaleW=') then
                  read(word,'(8x,f16.11)') scaleW
            else if(word(1:8).eq.'-scaleM=') then
                  read(word,'(8x,f16.11)') scaleM
            else if(word(1:8).eq.'-scaleS=') then
                  read(word,'(8x,f16.11)') scaleS
            end if
      end do

      ! REPORT LOG -----------------------------------------------------

      do idisp = 0,1
      write(idisp,*) ' '
      write(idisp,*) 'SCALE FOR WEIGHTS   (GIBBS)     : ', scaleW
      write(idisp,*) 'SCALE FOR MEANS     (GIBBS)     : ', scaleM
      write(idisp,*) 'SCALE FOR VARIANCES (GIBBS)     : ', scaleS
      write(idisp,*) '   '
      end do

! SCALE THE FIXED DIMENSION UPDATES (AIS) ==============================

      ! DEFAULTS
      scaleWais = scaleW
      scaleMais = scaleM
      scaleSais = scaleS

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:12).eq.'-ais-scaleW=') then
                  read(word,'(12x,f16.11)') scaleWais
            else if(word(1:12).eq.'-ais-scaleM=') then
                  read(word,'(12x,f16.11)') scaleMais
            else if(word(1:12).eq.'-ais-scaleS=') then
                  read(word,'(12x,f16.11)') scaleSais
            end if
      end do

      ! REPORT LOG -----------------------------------------------------

      do idisp = 0,1
      write(idisp,*) ' '
      write(idisp,*) 'SCALE FOR WEIGHTS   (AIS)       : ', scaleWais
      write(idisp,*) 'SCALE FOR MEANS     (AIS)       : ', scaleMais
      write(idisp,*) 'SCALE FOR VARIANCES (AIS)       : ', scaleSais
      write(idisp,*) ' '
      end do

! INITIALIZE THE ESTIMATIONS ===========================================

      sumK = 0
      sumKsq = 0

      ExpAccPrW = 0.d0
      ExpAccPrM = 0.d0
      ExpAccPrS = 0.d0

      nSplit = 0
      nSplitNadj = 0
      nMerge = 0

      ExpAccPrSplit = 0.d0
      ExpAccPrMerge = 0.d0

      ExpAccPrWais = 0.d0
      ExpAccPrMais = 0.d0
      ExpAccPrSais = 0.d0

! SET THE SEEDS ========================================================

      ! DEFAULT MODEL
      k = int(0.5d0*(kmin+kmax))

      ! ARGUMENTS
      do ia = 1,iargc()
            call getarg(ia,word)
            if(word(1:7).eq.'-kseed=') then
                  read(word,'(7x,i3)') k
            end if
      end do

      ! PARAMETERS
      beta = 0.02d0*dble(ymax-ymin)**2
      if (j.eq.1) then
            w(1) = 1.d0
            m(1) = 0.5d0*(ymin+ymax)
            s(1) = beta/(alpha-1)
      else
            do j = 1,k
                  w(j) = 1.d0/k
                  m(j) = ymin +(j-1.d0)/(k-1)*(ymax-ymin)
                  s(j) = beta/(alpha-1)
            end do
      end if

      ! REPORT LOG -----------------------------------------------------
      do idisp = 0,1
      write(idisp,*) ' '
      write(idisp,*) 'SEEDS'
      write(idisp,*) '====='
      write(idisp,*) 'NUMBER OF COMPONENTS', k
      write(idisp,*) '   WEIGHTS         MEANS           VARIANCE'
      do j = 1,k
            write(idisp,'(30(f15.10, x))') w(j), m(j), s(j)
      end do
      write(idisp,*) ' '
      write(idisp,*) '   BETA'
      write(idisp,'(f15.10, x)') beta
      write(idisp,*) ' '
      end do

! RUN THE RJMCMC SAMPLER ===============================================

      do iter = 1-nburnin,nsweeps

      ! PRINT THE ITERATION

            if (mod(dble(nsweeps-iter),dble(nsweeps/10)).eq.0) then
                  write(0,'(f5.1,X,$)') (nsweeps-iter)/(nsweeps/10.0d0)
            end if

      ! GIBBS BLOCK ====================================================
      ! ================================================================

            ! RECORD THE MODELS

            if (iter.gt.0) then
                  k1 = k
                  k2 = k
            end if

            ! UPDATE

            call FixedMoveSweep(k, &
                              w(1:k),m(1:k),s(1:k),beta,&
                              delta,xi,kappa,alpha,ge,ha,&
                              en,y(1:en),&
                              accPrW,accPrM,accPrS,&
                              scaleW,scaleM,scaleS)

            ! ESTIMATION OF THE ACCEPTANCE PROBABILITIES (FIXED)

            if (iter.gt.0) then
                  ExpAccPrW = ExpAccPrW+AccPrW
                  ExpAccPrM = ExpAccPrM+AccPrM
                  ExpAccPrS = ExpAccPrS+AccPrS
            end if

      ! Split/Merge BLOCK ==============================================
      ! ================================================================

      if (Tau.ge.1) then

      call rnguniform(u)

      if (qsplit(k).ge.u) then                              ! SPLIT MOVE

            ! RECORD THE MODELS
            if (iter.gt.0) then
                  k1 = k
                  k2 = k+1
            end if

      ! AIS DIMENSION MATCHING STEP TIME (t=0) =========================

      ! COMPUTE THE INVERSE TEMPERATURE

            t = 0
            gt = dble(t)/dble(Tau)
            Dgt = 1.d0/dble(Tau)

      ! CHOOSE THE COMPONENT TO BE SPLITED

            call rnguniform(u)
            j0 = 1+int((k)*u)
            wj0 = w(j0)
            mj0 = m(j0)
            sj0 = s(j0)

      ! GENERATE THE AUXILIARY VALUES

            call rngbeta(waux,2.d0,2.d0)
            call rngbeta(maux,2.d0,2.d0)
            call rngbeta(saux,1.d0,1.d0)

      ! COMPUTE THE DIMENSION MATCHING TRANSFORMATION
            wj1 = wj0*waux
            wj2 = wj0*(1-waux)
            mj1 = mj0 -maux*sqrt(sj0)*sqrt(wj2/wj1)
            mj2 = mj0 +maux*sqrt(sj0)*sqrt(wj1/wj2)
            ! C.I. -CHECK
            qadjacent = .true.
            if (j0.gt.1) then
                  if (m(j0-1).gt.mj1) then
                        qadjacent = .false.
                  end if
            end if
            if (j0.lt.k) then
                  if (mj2.gt.m(j0+1)) then
                        qadjacent = .false.
                  end if
            end if
            sj1 =     saux*(1-maux**2)*sj0*wj0/wj1
            sj2 = (1-saux)*(1-maux**2)*sj0*wj0/wj2

      ! LOAD THE AIS VARIABLES
            kais = k
            betaais = beta
            do j = 1,(j0-1)
                  wais(j) = w(j)
                  mais(j) = m(j)
                  sais(j) = s(j)
            end do
            wais((j0)) = wj1
            mais((j0)) = mj1
            sais((j0)) = sj1
            wais((j0+1)) = wj2
            mais((j0+1)) = mj2
            sais((j0+1)) = sj2
            do j = (j0+1)+1,k+1
                  wais(j) = w(j-1)
                  mais(j) = m(j-1)
                  sais(j) = s(j-1)
            end do

            ! COMPUTE THE LOGIW AT (t)

            if (qadjacent) then
                  call logImportanceWeightSplit(logIWt, &
                                          kais, &
                                          wais(1:kais+1),&
                                          mais(1:kais+1),&
                                          sais(1:kais+1), &
                                          j0, &
                                          betaais, &
                                          delta,xi,kappa,alpha, &
                                          en,y(1:en))
                  logAIW = Dgt*logIWt
            else
                  logAIW = -huge(logAIW)
            end if

      ! AIS STEP =======================================================

            if (qadjacent) then
            do t = 1,Tau-1

            ! COMPUTE THE INVERSE TEMPERATURE

                  gt = dble(t)/dble(Tau)
                  Dgt = 1.d0/dble(Tau)

            ! DO THE AIS SWEEP

                  call AisSplitSweep(kais, &
                                    wais(1:kais+1),delta, &
                                    mais(1:kais+1),xi,kappa, &
                                    sais(1:kais+1),alpha,&
                                    betaais,ge,ha,j0,en,y(1:en),gt,&
                                    AccPrWais,AccPrMais,AccPrSais,&
                                    scaleWais,scaleMais,scaleSais)

            ! COMPUTE THE LOGIW AT (t)

                  call logImportanceWeightSplit(logIWt, &
                                          kais, &
                                          wais(1:kais+1),&
                                          mais(1:kais+1),&
                                          sais(1:kais+1), &
                                          j0, &
                                          betaais, &
                                          delta,xi,kappa,alpha, &
                                          en,y(1:en))
                  logAIW = logAIW +Dgt*logIWt

            ! RECORD THE COUNTERS

                  if (iter.gt.0) then
                  ! ACCEPTANCE PROBABILITIES
                        ExpAccPrWais = ExpAccPrWais+AccPrWais
                        ExpAccPrMais = ExpAccPrMais+AccPrMais
                        ExpAccPrSais = ExpAccPrSais+AccPrSais
                  end if

            end do
            end if

      ! AIS ACCEPT / REJECT STEP =======================================

            if (qadjacent) then
                  AccPrRJ = min(1.d0,qmerge(k+1)/qsplit(k)*exp(logAIW))
            else
                  AccPrRJ = 0.d0
            end if

            call rnguniform(u)
            if (AccPrRJ.gt.u) then
                  k = k+1
                  beta = betaais
                  do j = 1,k
                        w(j) = wais(j)
                        m(j) = mais(j)
                        s(j) = sais(j)
                  end do
            end if

      ! RECORD THE MOVE ================================================

            ! COUNTERS
            if (iter.gt.0) then
                  nSplit = nSplit+1
                  if (.not.qadjacent) nSplitNadj = nSplitNadj+1
                  ExpAccPrSplit = ExpAccPrSplit + AccPrRJ
            end if

      elseif (qsplit(k)+qmerge(k).ge.u) then                ! MERGE MOVE

            ! RECORD THE MODELS
            if (iter.gt.0) then
                  k1 = k
                  k2 = k-1
            end if

      ! AIS DIMENSION MATCHING STEP TIME (t=Tau) =======================

            ! COMPUTE THE INVERSE TEMPERATURE

            t = Tau
            gt = dble(t)/dble(Tau)
            Dgt = 1.d0/dble(Tau)

            ! CHOOSE TO MERGE J0 -TH AND J0+1 -TH COMPONENTS

            call rnguniform(u)
            j0 = 1+int((k-1)*u)

            ! LOAD THE AIS VARIABLES

            kais = k-1
            betaais = beta
            do j = 1,k
                  wais(j) = w(j)
                  mais(j) = m(j)
                  sais(j) = s(j)
            end do

            ! COMPUTE THE LOGIW AT (t=0)

            call logImportanceWeightSplit(logIWt, &
                                          kais, &
                                          wais(1:kais+1),&
                                          mais(1:kais+1),&
                                          sais(1:kais+1), &
                                          j0, &
                                          betaais, &
                                          delta,xi,kappa,alpha, &
                                          en,y(1:en))
            logAIW = Dgt*(-logIWt)

      ! AIS STEP =======================================================

            do t = Tau-1,1,-1

            ! COMPUTE THE INVERSE TEMPERATURE AT (t)

                  gt = dble(t)/dble(Tau)
                  Dgt = 1.d0/dble(Tau)

            ! DO THE AIS SWEEP AT (t)

                  call AisSplitSweep(kais, &
                                    wais(1:kais+1),delta, &
                                    mais(1:kais+1),xi,kappa, &
                                    sais(1:kais+1),alpha,&
                                    betaais,ge,ha,j0,en,y(1:en),gt,&
                                    AccPrWais,AccPrMais,AccPrSais,&
                                    scaleWais,scaleMais,scaleSais)

            ! COMPUTE THE LOGIW AT (t)

                  call logImportanceWeightSplit(logIWt, &
                                                kais, &
                                                wais(1:kais+1),&
                                                mais(1:kais+1),&
                                                sais(1:kais+1), &
                                                j0, &
                                                betaais, &
                                                delta,xi,kappa,alpha, &
                                                en,y(1:en))
                  logAIW = logAIW +Dgt*(-logIWt)

            ! RECORD THE COUNTERS

                  if (iter.gt.0) then
                  ! ACCEPTANCE PROBABILITIES
                        ExpAccPrWais = ExpAccPrWais+AccPrWais
                        ExpAccPrMais = ExpAccPrMais+AccPrMais
                        ExpAccPrSais = ExpAccPrSais+AccPrSais
                  end if

            end do

      ! AIS ACCEPT / REJECT STEP =======================================

            AccPrRJ = min(1.d0,qsplit(k-1)/qmerge(k)*exp(logAIW) )

            call rnguniform(u)
            if (AccPrRJ.gt.u) then
                  k = k-1
                  beta = betaais
                  do j = 1,j0-1
                        w(j) = wais(j)
                        m(j) = mais(j)
                        s(j) = sais(j)
                  end do
                  w(j0) = wais((j0))+wais((j0+1))
                  m(j0) = (wais((j0))*mais((j0))+wais((j0+1))*mais((j0+1)))/w(j0)
                  s(j0) = ( &
                              wais((j0))*(mais((j0))**2+sais((j0))) &
                              +wais((j0+1))*(mais((j0+1))**2+sais((j0+1))) &
                              ) &
                              /w(j0)-m(j0)**2
                  do j = j0+1,k
                        w(j) = wais(j+1)
                        m(j) = mais(j+1)
                        s(j) = sais(j+1)
                  end do
            end if

      ! RECORD THE MOVE ================================================

            ! COUNTERS
            if (iter.gt.0) then
                  nMerge = nMerge+1
                  ExpAccPrMerge = ExpAccPrMerge + AccPrRJ
            end if

      end if

      end if ! [if (Tau.ge.1)] ! LOCAL Split & Merge RJ MOVE

      ! SAVE ===========================================================

      if (iter.gt.0) then

      ! STORE
            write(3,'(i2)') k
            write(4,'(100(f20.11, x))') (w(j),j=1,k)
            write(5,'(100(f20.11, x))') (m(j),j=1,k)
            write(6,'(100(f20.11, x))') (s(j),j=1,k)
            write(7,'(f20.11)') beta

            write(8,'(i3, x, i3)') k1, k2

            write(9,*) logAIW

            write(10,*) AccPrRJ

      ! ESTIMATION OF mean(k)
            sumK = sumK+k
            sumKsq = sumKsq+k**2

      end if

      end do ! THE SWEEP

! REPORT LOG ===========================================================

      do idisp = 0,1
      write(idisp,*) '   '
      write(idisp,*) '   '
      write(idisp,'( "Posterior k mean",f15.10,"(",f15.10,")" )') &
             sumK/dble(nsweeps),sqrt(sumKsq/dble(nsweeps)-(sumK/dble(nsweeps))**2)

      write(idisp,*) '   '
      write(idisp,'("Expected acceptance probabilities for the fixed moves")')
      write(idisp,'("Expected acceptance probability (weights move)        : ",f15.10)') &
             ExpAccPrW/(nsweeps)
      write(idisp,'("Expected acceptance probability (means move)          : ",f15.10)') &
             ExpAccPrM/(nsweeps)
      write(idisp,'("Expected acceptance probability (variances move)      : ",f15.10)') &
             ExpAccPrS/(nsweeps)
      end do


      do idisp = 0,1
      if (Tau.gt.1) then
      write(idisp,*) '   '
      write(idisp,'("Expected acceptance probabilities for the fixed moves of AIS")')
      write(idisp,'("Expected acceptance probability (weights move)        : ",f15.10)') &
             ExpAccPrWais/((nSplit-nSplitNadj+nMerge)*(Tau-1))
      write(idisp,'("Expected acceptance probability (means move)          : ",f15.10)') &
             ExpAccPrMais/((nSplit-nSplitNadj+nMerge)*(Tau-1))
      write(idisp,'("Expected acceptance probability (variances move)      : ",f15.10)') &
             ExpAccPrSais/((nSplit-nSplitNadj+nMerge)*(Tau-1))
      end if
      end do

      do idisp = 0,1
      if (Tau.ge.1) then
      write(idisp,*) '   '
      write(idisp,'("Expected acceptance probabilities for the RJ moves")')
      write(idisp,'("Expected acceptance probability (Split move)          : ",f15.10)') &
             ExpAccPrSplit/(nSplit)
      write(idisp,'("Expected acceptance probability (Merge move)          : ",f15.10)') &
             ExpAccPrMerge/(nMerge)
      write(idisp,'("Expected acceptance probability (Split/Merge move)    : ",f15.10)') &
             (ExpAccPrSplit+ExpAccPrMerge)/(nSplit+nMerge)
      end if
      end do

      do idisp = 0,1
      if (Tau.ge.1) then
      write(idisp,*) '   '
      write(idisp,'("Rejections at the split jump")')
      write(idisp,'("Expected percentage of rejections                     : ",f15.10)') &
             dble(nSplitNadj)/nSplit
      write(idisp,'("Number of rejections                                  : ",i15)') &
             nSplitNadj
      write(idisp,'("Number of split moves                                 : ",i15)') &
             nSplit
      end if
      end do

! ======================================================================

end program aisrjnmix_pro


