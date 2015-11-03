! 2010/09/01

!      INFORMATIONS:
!      =============
!      THE SPLIT-MERGE MOVE
!      PRIORS                 : DU(INDEX)
!      RELABELLING            : IC
!      GENERATE/KILL          : RICHARDSON AND GREEN 1997
!      MARGINAL MODEL
!      =================================================================

module splitmergemove_mod

      implicit none

      private
      public      :: logImportanceWeightSplit
      public      :: AisSplitSweep

contains

!=======================================================================
! ACCEPTANCE RATIO
!=======================================================================

      subroutine logImportanceWeightSplit(logAIW, &
                                          k,w,m,s,j0,&
                                          beta,delta,xi,kappa,alpha,&
                                          en,y)

      ! LOG FORM

            implicit none

            ! DECLARE THE ARGUMENTS

            double precision, intent(out) :: logAIW
            integer, intent(in)           :: k
            integer, intent(in)           :: en
            double precision, intent(in)  :: w(k+1)
            double precision, intent(in)  :: m(k+1)
            double precision, intent(in)  :: s(k+1)
            integer, intent(in)           :: j0
            double precision, intent(in)  :: beta
            double precision, intent(in)  :: delta
            double precision, intent(in)  :: alpha
            double precision, intent(in)  :: xi
            double precision, intent(in)  :: kappa
            double precision, intent(in)  :: y(en)

            ! DECLARE THE PARAMETERS
            double precision  :: pi
            parameter (pi=3.14159265358979d0)

            ! DECLARE THE LOCAL FUNCTIONS
            integer           :: i,j
            double precision  :: yi
            double precision  :: wj0,mj0,sj0
            double precision  :: wj1,wj2,mj1
            double precision  :: mj2,sj1,sj2
            double precision  :: waux,maux,saux
            double precision  :: CJc,CJ0,CJ1,CJ2
            double precision  :: dlgama

            wj1 = w((j0))
            mj1 = m((j0))
            sj1 = s((j0))
            wj2 = w((j0+1))
            mj2 = m((j0+1))
            sj2 = s((j0+1))

            wj0 = wj1+wj2
            mj0 = (wj1*mj1+wj2*mj2)/wj0
            sj0 = (wj1*(mj1**2+sj1)+wj2*(mj2**2+sj2))/wj0-mj0**2

            waux = wj1/wj0
            maux = (mj2-mj0)*sqrt(wj2/wj1)/sqrt(sj0)
            saux = sj1*wj1/((1-maux**2)*sj0*wj0)

            logAIW = 0.d0

            do i = 1,en
                  yi = y(i)
                  CJc = 0.d0
                  do j = 1,(j0)-1
                        CJc = CJc &
                              +w(j)/sqrt(s(j))*exp(-0.5d0*(yi-m(j))**2/s(j))
                  end do
                  CJ0 = wj0/sqrt(sj0)*exp(-0.5d0*(yi-mj0)**2/sj0)
                  CJ1 = wj1/sqrt(sj1)*exp(-0.5d0*(yi-mj1)**2/sj1)
                  CJ2 = wj2/sqrt(sj2)*exp(-0.5d0*(yi-mj2)**2/sj2)
                  do j = (j0+1)+1,k+1
                        CJc = CJc &
                              +w(j)/sqrt(s(j))*exp(-0.5d0*(yi-m(j))**2/s(j))
                  end do
                  logAIW = logAIW &
                              +log(CJc+CJ1+CJ2)-log(CJc+CJ0)
            end do
            
            logAIW = logAIW &
                  ! weights
                        +dlgama(delta*k+delta)-dlgama(delta*k)-dlgama(delta) &
                        +(delta-1)*log(wj1) &
                        +(delta-1)*log(wj2) &
                        -(delta-1)*log(wj0) &
                  ! mean
                        +0.5d0*(log(kappa)-log(2*pi)) &
                        -kappa/2*((mj1-xi)**2+(mj2-xi)**2-(mj0-xi)**2) &
                        +log(k+1.d0) &
                  ! variance
                        +alpha*log(beta) &
                        -dlgama(alpha) &
                        -(alpha+1)*(log(sj1)+log(sj2)-log(sj0)) &
                        -beta*(1/sj1+1/sj2-1/sj0) &
                  ! auxiliaries
                        -log(36.d0)-log(waux)-log(1-waux)-log(maux)-log(1-maux) &
                  ! Jucobian
                        +log(wj0)+log(abs(mj1-mj2))+log(sj1)+log(sj2) &
                        -log(maux)-log(1-maux**2)-log(saux)-log(1-saux)-log(sj0)

      end subroutine logImportanceWeightSplit

!=======================================================================
! AIS SPLIT SWEEP
!=======================================================================

      subroutine AisSplitSweep(k,w,delta,&
                                    m,xi,kappa,&
                                    s,alpha,&
                                    beta,ge,ha,&
                                    j0,en,y,gt,&
                                    AccPrW,AccPrM,AccPrS,&
                                    sclw,sclm,scls)

            implicit none

            ! DECLARE THE ARGUMENTS
            integer, intent(in)                 :: k
            double precision, intent(inout)     :: w(k+1)
            double precision, intent(in)        :: delta
            double precision, intent(inout)     :: m(k+1)
            double precision, intent(in)        :: xi,kappa
            double precision, intent(inout)     :: s(k+1)
            double precision, intent(in)        :: alpha
            double precision, intent(inout)     :: beta
            double precision, intent(in)        :: ge,ha
            integer, intent(inout)              :: j0
            integer, intent(in)                 :: en
            double precision, intent(in)        :: y(en)
            double precision, intent(in)        :: gt
            double precision, intent(out)       :: accPrW
            double precision, intent(out)       :: accPrM
            double precision, intent(out)       :: accPrS
            double precision, intent(in)        :: sclw
            double precision, intent(in)        :: sclm
            double precision, intent(in)        :: scls

            ! DECLARE THE LOCAL VARIABLES
            integer           :: j,i
            double precision  :: yi
            double precision  :: dvar,dvar0,dvar1,dvar2
            double precision  :: dvarN,dvarN0,dvarN1,dvarN2
            double precision  :: vec(k+1)
            double precision  :: wj0,mj0,sj0
            double precision  :: wj1,mj1,sj1
            double precision  :: wj2,mj2,sj2
            double precision  :: aux1,aux2,aux3
            double precision  :: wNj0,mNj0,sNj0
            double precision  :: wNj1,mNj1,sNj1
            double precision  :: wNj2,mNj2,sNj2
            double precision  :: auxN1,auxN2,auxN3
            double precision  :: logRho,logRhoN
            double precision  :: logPropN,logProp
            double precision  :: u
            logical           :: isOrdered

            double precision  :: Snumer
            double precision  :: pr(k)


            integer           :: nbl
            parameter         ( nbl = 5 )
            integer           :: ibl,rbl(nbl)

            ! GENERATE A RANDOM ORDER FOR THE BLOCKS
            do j = 1,nbl
                  rbl(j) = j
            end do
            call genprm(rbl(1:nbl),nbl)

            ! RANDOM PERMUTATION SCAN MCMC SWEEP
      do ibl = 1,nbl

      if (rbl(ibl).eq.1) then
      ! WEIGHT BLOCK ; FORWARD
      ! ======================

            ! COMPUTE THE PROPOSALS

            dvar1 = w(k+1)
            do j = 1, k
                  call rngnormal(vec(j:j),1)
                  vec(j) = w(j)*exp(sclw*vec(j))
                  dvar1 = dvar1+vec(j)
            end do
            dvar2 = 0.d0
            do j = 1,k
                  vec(j) = vec(j)/dvar1
                  dvar2 = dvar2 +vec(j)
            end do
            vec(k+1) = 1-dvar2

            ! COMPUTE THE ACCEPTANCE PROBABILITY

            wj1 = w((j0))
            mj1 = m((j0))
            sj1 = s((j0))
            wj2 = w((j0+1))
            mj2 = m((j0+1))
            sj2 = s((j0+1))
            wj0 = wj1+wj2
            mj0 = (wj1*mj1+wj2*mj2)/wj0
            sj0 = (wj1*(mj1**2+sj1)+wj2*(mj2**2+sj2))/wj0-mj0**2
            aux1 = wj1/wj0
            aux2 = (mj2-mj0)*sqrt(wj2/wj1)/sqrt(sj0)
            aux3 = sj1*wj1/((1-aux2**2)*sj0*wj0)

            wNj1 = vec((j0))
            wNj2 = vec((j0+1))
            wNj0 = wNj1+wNj2
            mNj0 = (wNj1*mj1+wNj2*mj2)/wNj0
            sNj0 = (wNj1*(mj1**2+sj1)+wNj2*(mj2**2+sj2))/wNj0-mNj0**2
            auxN1 = wNj1/wNj0
            auxN2 = (mj2-mNj0)*sqrt(wNj2/wNj1)/sqrt(sNj0)
            auxN3 = sj1*wNj1/((1-auxN2**2)*sNj0*wNj0)

            logRho = 0.d0
            logRhoN = 0.d0

            ! LIKELIHOOD
            do i = 1,en
                  yi = y(i)
                  dvar = 0.d0
                  dvarN = 0.d0
                  do j = 1,(j0)-1
                        dvar = dvar +w(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/s(j))
                        dvarN = dvarN +vec(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/s(j))
                  end do
                  dvar0 = wj0/sqrt(sj0) &
                              *exp(-0.5d0*(yi-mj0)**2/sj0)
                  dvarN0 = wNj0/sqrt(sNj0) &
                              *exp(-0.5d0*(yi-mNj0)**2/sNj0)
                  dvar1 = wj1/sqrt(sj1) &
                              *exp(-0.5d0*(yi-mj1)**2/sj1)
                  dvarN1 = wNj1/sqrt(sj1) &
                              *exp(-0.5d0*(yi-mj1)**2/sj1)
                  dvar2 = wj2/sqrt(sj2) &
                              *exp(-0.5d0*(yi-mj2)**2/sj2)
                  dvarN2 = wNj2/sqrt(sj2) &
                              *exp(-0.5d0*(yi-mj2)**2/sj2)
                  do j = (j0+1)+1,k+1
                        dvar = dvar +w(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/s(j))
                        dvarN = dvarN +vec(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/s(j))
                  end do
                  logRho = logRho &
                              +(1-gt)*log(dvar+dvar0) &
                                    +gt*log(dvar+dvar1+dvar2)
                  logRhoN = logRhoN &
                              +(1-gt)*log(dvarN+dvarN0) &
                                    +gt*log(dvarN+dvarN1+dvarN2)
            end do
            ! WEIGHTS PRIOR
            do j = 1,(j0)-1
                  logRho = logRho +(delta-1)*log(w(j))
                  logRhoN = logRhoN +(delta-1)*log(vec(j))
            end do
            logRho = logRho +(1-gt)*(delta-1)*log(wj0) &
                              +gt*(delta-1)*log(wj1) &
                              +gt*(delta-1)*log(wj2)
            logRhoN = logRhoN +(1-gt)*(delta-1)*log(wNj0) &
                              +gt*(delta-1)*log(wNj1) &
                              +gt*(delta-1)*log(wNj2)
            do j = (j0+1)+1,k+1
                  logRho = logRho +(delta-1)*log(w(j))
                  logRhoN = logRhoN +(delta-1)*log(vec(j))
            end do
            ! MEANS PRIOR
            logRho = logRho -(1-gt)*0.5d0*(mj0-xi)**2*kappa
            logRhoN = logRhoN -(1-gt)*0.5d0*(mNj0-xi)**2*kappa
            ! VARIANCES PRIOR
            logRho = logRho -(1-gt)*(alpha+1)*log(sj0)-(1-gt)*beta/sj0
            logRhoN = logRhoN -(1-gt)*(alpha+1)*log(sNj0)-(1-gt)*beta/sNj0
            ! AUXILIARY VALUES
            logRho = logRho &
                        +(1-gt)*log(aux1)+(1-gt)*log(1-aux1) &
                        +(1-gt)*log(aux2)+(1-gt)*log(1-aux2)
            logRhoN = logRhoN &
                        +(1-gt)*log(auxN1)+(1-gt)*log(1-auxN1) &
                        +(1-gt)*log(auxN2)+(1-gt)*log(1-auxN2)
            ! JACOBIAN
            logRho = logRho &
                        -(1-gt)*( &
                              log(wj0) &
                              -log(aux2)-log(1-aux2**2) &
                              -log(aux3)-log(1-aux3) &
                              -log(sj0) &
                              )
            logRhoN = logRhoN &
                        -(1-gt)*( &
                              log(wNj0) &
                              -log(auxN2)-log(1-auxN2**2) &
                              -log(auxN3)-log(1-auxN3) &
                              -log(sNj0) &
                              )

            ! PROPOSALS
            logProp = 0.d0
            logPropN = 0.d0
            do j = 1,k+1
                  logProp = logProp -log(w(j))
                  logPropN = logPropN -log(vec(j))
            end do

            ! ACCEPT/REJECT
            accPrW = min(1.d0, &
                        exp(logRhoN-logRho+logProp-logPropN) &
                        )
            call rnguniform(u)
            if (accPrW.gt.u) then
                  do j = 1,k+1
                        w(j) = vec(j)
                  end do
            end if

      elseif (rbl(ibl).eq.2) then
      ! MEAN BLOCK ; FORWARD
      ! ====================

            ! GENERATE PROPOSALS

            do j = 1,k+1
                  call rngnormal(vec(j:j),1)
                  vec(j) = m(j)+vec(j)*sclm
            end do

            ! CHECK THE ORDER OF THE PROPOSAL

            isOrdered = .true.
            do j = 1,k
                  if (vec(j).gt.vec(j+1)) then
                        isOrdered = .false.
                        exit
                  end if
            end do

            ! COMPUTE THE ACCEPTANCE PROBABILITY

            if (isOrdered) then

                  wj1 = w((j0))
                  mj1 = m((j0))
                  sj1 = s((j0))
                  wj2 = w((j0+1))
                  mj2 = m((j0+1))
                  sj2 = s((j0+1))
                  wj0 = wj1+wj2
                  mj0 = (wj1*mj1+wj2*mj2)/wj0
                  sj0 = (wj1*(mj1**2+sj1)+wj2*(mj2**2+sj2))/wj0-mj0**2
                  aux1 = wj1/wj0
                  aux2 = (mj2-mj0)*sqrt(wj2/wj1)/sqrt(sj0)
                  aux3 = sj1*wj1/((1-aux2**2)*sj0*wj0)

                  mNj1 = vec((j0))
                  mNj2 = vec((j0+1))
                  wNj0 = wj1+wj2
                  mNj0 = (wj1*mNj1+wj2*mNj2)/wNj0
                  sNj0 = (wj1*(mNj1**2+sj1)+wj2*(mNj2**2+sj2))/wNj0-mNj0**2
                  auxN1 = wj1/wNj0
                  auxN2 = (mNj2-mNj0)*sqrt(wj2/wj1)/sqrt(sNj0)
                  auxN3 = sj1*wj1/((1-auxN2**2)*sNj0*wNj0)

                  logRho = 0.d0
                  logRhoN = 0.d0

                  ! LIKELIHOOD
                  do i = 1,en
                        yi = y(i)
                        dvar = 0.d0
                        dvarN = 0.d0
                        do j = 1,(j0)-1
                              dvar = dvar +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(yi-m(j))**2/s(j))
                              dvarN = dvarN +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(yi-vec(j))**2/s(j))
                        end do
                        dvar0 = wj0/sqrt(sj0) &
                                    *exp(-0.5d0*(yi-mj0)**2/sj0)
                        dvarN0 = wNj0/sqrt(sNj0) &
                                    *exp(-0.5d0*(yi-mNj0)**2/sNj0)
                        dvar1 = wj1/sqrt(sj1) &
                                    *exp(-0.5d0*(yi-mj1)**2/sj1)
                        dvarN1 = wj1/sqrt(sj1) &
                                    *exp(-0.5d0*(yi-mNj1)**2/sj1)
                        dvar2 = wj2/sqrt(sj2) &
                                    *exp(-0.5d0*(yi-mj2)**2/sj2)
                        dvarN2 = wj2/sqrt(sj2) &
                                    *exp(-0.5d0*(yi-mNj2)**2/sj2)
                        do j = (j0+1)+1,k+1
                              dvar = dvar +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(yi-m(j))**2/s(j))
                              dvarN = dvarN +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(yi-vec(j))**2/s(j))
                        end do
                        logRho = logRho &
                                    +(1-gt)*log(dvar+dvar0) &
                                          +gt*log(dvar+dvar1+dvar2)
                        logRhoN = logRhoN &
                                    +(1-gt)*log(dvarN+dvarN0) &
                                          +gt*log(dvarN+dvarN1+dvarN2)
                  end do

                  ! WEIGHTS PRIOR
                  logRho = logRho +(1-gt)*(delta-1)*log(wj0)
                  logRhoN = logRhoN +(1-gt)*(delta-1)*log(wNj0)

                  ! MEANS PRIOR
                  do j = 1,(j0)-1
                        logRho = logRho -0.5d0*(m(j)-xi)**2*kappa
                        logRhoN = logRhoN -0.5d0*(vec(j)-xi)**2*kappa
                  end do
                  logRho = logRho &
                              -(1-gt)*0.5d0*(mj0-xi)**2*kappa &
                              -gt*0.5d0*(mj1-xi)**2*kappa &
                              -gt*0.5d0*(mj2-xi)**2*kappa
                  logRhoN = logRhoN &
                              -(1-gt)*0.5d0*(mNj0-xi)**2*kappa &
                              -gt*0.5d0*(mNj1-xi)**2*kappa &
                              -gt*0.5d0*(mNj2-xi)**2*kappa
                  do j = (j0+1)+1,k+1
                        logRho = logRho -0.5d0*(m(j)-xi)**2*kappa
                        logRhoN = logRhoN -0.5d0*(vec(j)-xi)**2*kappa
                  end do

                  ! VARIANCES PRIORS
                  logRho = logRho &
                              -(1-gt)*(alpha+1)*log(sj0) &
                              -(1-gt)*beta/sj0
                  logRhoN = logRhoN &
                              -(1-gt)*(alpha+1)*log(sNj0) &
                              -(1-gt)*beta/sNj0
                  ! AUXILIARY VALUES
                  logRho = logRho &
                              +(1-gt)*log(aux1) +(1-gt)*log(1-aux1) &
                              +(1-gt)*log(aux2) +(1-gt)*log(1-aux2)
                  logRhoN = logRhoN &
                              +(1-gt)*log(auxN1) +(1-gt)*log(1-auxN1) &
                              +(1-gt)*log(auxN2) +(1-gt)*log(1-auxN2)

                  ! JACOBIAN
                  logRho = logRho &
                              -(1-gt)*( &
                                    log(wj0)+log(abs(mj2-mj1)) &
                                    -log(aux2)-log(1-aux2**2) &
                                    -log(aux3)-log(1-aux3) &
                                    -log(sj0) &
                                    )
                  logRhoN = logRhoN &
                              -(1-gt)*( &
                                    log(wNj0)+log(abs(mNj2-mNj1)) &
                                    -log(auxN2)-log(1-auxN2**2) &
                                    -log(auxN3)-log(1-auxN3) &
                                    -log(sNj0) &
                                    )

                  ! COMPUTE PROPOSALS DENSITIES (up to the other proposal)
                  logProp = 0.d0
                  logPropN = 0.d0

                  ! ACCEPTANCE PROBABILITY
                  accPrM = min(1.d0, &
                              exp(logRhoN-logRho+logProp-logPropN) &
                              )

            else

                  accPrM = 0.d0

            end if

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (accPrM.gt.u) then
                  do j = 1,k+1
                        m(j) = vec(j)
                  end do
            end if

      elseif (rbl(ibl).eq.3) then
      ! VARIANCES BLOCK ; FORWARD
      ! =========================

            ! GENERATE PROPOSALS

            do j = 1,k+1
                  call rngnormal(vec(j:j),1)
                  vec(j) = s(j)*exp(vec(j)*scls)
            end do

            ! COMPUTE THE ACCEPTANCE PROBABILITY

            wj1 = w((j0))
            mj1 = m((j0))
            sj1 = s((j0))
            wj2 = w((j0+1))
            mj2 = m((j0+1))
            sj2 = s((j0+1))
            wj0 = wj1+wj2
            mj0 = (wj1*mj1+wj2*mj2)/wj0
            sj0 = (wj1*(mj1**2+sj1)+wj2*(mj2**2+sj2))/wj0-mj0**2
            aux1 = wj1/wj0
            aux2 = (mj2-mj0)*sqrt(wj2/wj1)/sqrt(sj0)
            aux3 = sj1*wj1/((1-aux2**2)*sj0*wj0)

            sNj1 = vec((j0))
            sNj2 = vec((j0+1))
            wNj0 = wj1+wj2
            mNj0 = (wj1*mj1+wj2*mj2)/wNj0
            sNj0 = (wj1*(mj1**2+sNj1)+wj2*(mj2**2+sNj2))/wNj0-mNj0**2
            auxN1 = wj1/wNj0
            auxN2 = (mj2-mNj0)*sqrt(wj2/wj1)/sqrt(sNj0)
            auxN3 = sNj1*wj1/((1-auxN2**2)*sNj0*wNj0)

            logRho = 0.d0
            logRhoN = 0.d0

            ! LIKELIHOOD
            do i = 1,en
                  yi = y(i)
                  dvar = 0.d0
                  dvarN = 0.d0
                  do j = 1,(j0)-1
                        dvar = dvar +w(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/s(j))
                        dvarN = dvarN +w(j)/sqrt(vec(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/vec(j))
                  end do
                  dvar0 = wj0/sqrt(sj0) &
                              *exp(-0.5d0*(yi-mj0)**2/sj0)
                  dvarN0 = wNj0/sqrt(sNj0) &
                              *exp(-0.5d0*(yi-mNj0)**2/sNj0)
                  dvar1 = wj1/sqrt(sj1) &
                              *exp(-0.5d0*(yi-mj1)**2/sj1)
                  dvarN1 = wj1/sqrt(sNj1) &
                              *exp(-0.5d0*(yi-mj1)**2/sNj1)
                  dvar2 = wj2/sqrt(sj2) &
                              *exp(-0.5d0*(yi-mj2)**2/sj2)
                  dvarN2 = wj2/sqrt(sNj2) &
                              *exp(-0.5d0*(yi-mj2)**2/sNj2)
                  do j = (j0+1)+1,k+1
                        dvar = dvar +w(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/s(j))
                        dvarN = dvarN +w(j)/sqrt(vec(j)) &
                                    *exp(-0.5d0*(yi-m(j))**2/vec(j))
                  end do
                  logRho = logRho &
                              +(1-gt)*log(dvar+dvar0) &
                                    +gt*log(dvar+dvar1+dvar2)
                  logRhoN = logRhoN &
                              +(1-gt)*log(dvarN+dvarN0) &
                                    +gt*log(dvarN+dvarN1+dvarN2)
            end do
            ! WEIGHTS PRIOR
            logRho = logRho +(1-gt)*(delta-1)*log(wj0)
            logRhoN = logRhoN +(1-gt)*(delta-1)*log(wNj0)
            ! MEANS PRIOR
            logRho = logRho &
                        -(1-gt)*0.5d0*(mj0-xi)**2*kappa 
            logRhoN = logRhoN &
                        -(1-gt)*0.5d0*(mNj0-xi)**2*kappa 
            ! VARIANCES PRIOR
            do j = 1,(j0)-1
                  logRho = logRho &
                              -(alpha+1)*log(s(j)) -beta/s(j)
                  logRhoN = logRhoN &
                              -(alpha+1)*log(vec(j)) -beta/vec(j)
            end do
            logRho = logRho &
                        -(1-gt)*(alpha+1)*log(sj0) -(1-gt)*beta/sj0 &
                        -gt*(alpha+1)*log(sj1) -gt*beta/sj1 &
                        -gt*(alpha+1)*log(sj2) -gt*beta/sj2
            logRhoN = logRhoN &
                        -(1-gt)*(alpha+1)*log(sNj0) -(1-gt)*beta/sNj0 &
                        -gt*(alpha+1)*log(sNj1) -gt*beta/sNj1 &
                        -gt*(alpha+1)*log(sNj2) -gt*beta/sNj2
            do j = (j0+1)+1,k+1
                  logRho = logRho &
                              -(alpha+1)*log(s(j)) -beta/s(j)
                  logRhoN = logRhoN &
                              -(alpha+1)*log(vec(j)) -beta/vec(j)
            end do
            ! AUXILIARY VALUES
            logRho = logRho &
                        +(1-gt)*log(aux1) +(1-gt)*log(1-aux1) &
                        +(1-gt)*log(aux2) +(1-gt)*log(1-aux2)
            logRhoN = logRhoN &
                        +(1-gt)*log(auxN1) +(1-gt)*log(1-auxN1) &
                        +(1-gt)*log(auxN2) +(1-gt)*log(1-auxN2)
            ! JACOBIAN
            logRho = logRho &
                        -(1-gt)*( &
                              log(wj0)+log(sj1)+log(sj2) &
                              -log(aux2)-log(1-aux2**2) &
                              -log(aux3)-log(1-aux3) &
                              -log(sj0) &
                              )
            logRhoN = logRhoN &
                        -(1-gt)*( &
                              log(wNj0)+log(sNj1)+log(sNj2) &
                              -log(auxN2)-log(1-auxN2**2) &
                              -log(auxN3)-log(1-auxN3) &
                              -log(sNj0) &
                              )

            ! COMPUTE PROPOSALS DENSITIES (up to the other proposal)
            logProp = 0.d0
            logPropN = 0.d0
            do j = 1,k+1
                  logProp = logProp -log(s(j))
                  logPropN = logPropN -log(vec(j))
            end do

            ! COMPUTE THE ACCEPTANCE PROBABILITY
            accPrS = min(1.d0, &
                        exp(logRhoN-logRho+logProp-logPropN) &
                        )

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (accPrS.gt.u) then
                  do j = 1,k+1
                        s(j) = vec(j)
                  end do
            end if

      elseif (rbl(ibl).eq.4) then
      ! AIS BETA BLOCK ; FORWARD
      ! ========================

            wj1 = w((j0))
            mj1 = m((j0))
            sj1 = s((j0))
            wj2 = w((j0+1))
            mj2 = m((j0+1))
            sj2 = s((j0+1))
            wj0 = wj1+wj2
            mj0 = (wj1*mj1+wj2*mj2)/wj0
            sj0 = (wj1*(mj1**2+sj1)+wj2*(mj2**2+sj2))/wj0-mj0**2

            dvar1 = ge +(k+gt)*alpha
            dvar2 = ha
            do j = 1,(j0)-1
                  dvar2 = dvar2 + 1.0d0/s(j)
            end do
            dvar2 = dvar2 +(1-gt)/sj0 +gt/sj1 +gt/sj2
            do j = (j0+1)+1,k+1
                  dvar2 = dvar2 +1.0d0/s(j)
            end do
            call rnggamma(beta,dvar1)
            beta = beta/dvar2

      elseif (rbl(ibl).eq.5) then
      ! AIS J0 BLOCK ; FORWARD
      ! ======================

            Snumer = 0.d0

            do j0 = 1,k

                  ! COMPUTE THE EXEPTIONS

                  wj1 = w((j0))
                  mj1 = m((j0))
                  sj1 = s((j0))
                  wj2 = w((j0+1))
                  mj2 = m((j0+1))
                  sj2 = s((j0+1))
                  wj0 = wj1+wj2
                  mj0 = (wj1*mj1+wj2*mj2)/wj0
                  sj0 = (wj1*(mj1**2+sj1)+wj2*(mj2**2+sj2))/wj0-mj0**2
                  aux1 = wj1/wj0
                  aux2 = (mj2-mj0)*sqrt(wj2/wj1)/sqrt(sj0)
                  aux3 = sj1*wj1/((1-aux2**2)*sj0*wj0)

                  ! COMPUTE THE NUMERATOR OF THE PROBABILITY

                  logRho = 0.d0

                  ! LIKELIHOOD
                  do i = 1,en
                        yi = y(i)
                        dvar = 0.d0
                        do j = 1,(j0)-1
                              dvar = dvar +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(yi-m(j))**2/s(j))
                        end do
                        dvar0 = wj0/sqrt(sj0) &
                                    *exp(-0.5d0*(yi-mj0)**2/sj0)
                        dvar1 = wj1/sqrt(sj1) &
                                    *exp(-0.5d0*(yi-mj1)**2/sj1)
                        dvar2 = wj2/sqrt(sj2) &
                                    *exp(-0.5d0*(yi-mj2)**2/sj2)
                        do j = (j0+1)+1,k+1
                              dvar = dvar +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(yi-m(j))**2/s(j))
                        end do
                        logRho = logRho &
                                    +(1-gt)*log(dvar+dvar0) &
                                          +gt*log(dvar+dvar1+dvar2)
                  end do

                  ! WEIGHTS PRIOR
                  do j = 1,(j0)-1
                        logRho = logRho +(delta-1)*log(w(j))
                  end do
                  logRho = logRho +(1-gt)*(delta-1)*log(wj0) &
                                    +gt*(delta-1)*log(wj1) &
                                    +gt*(delta-1)*log(wj2)
                  do j = (j0+1)+1,k+1
                        logRho = logRho +(delta-1)*log(w(j))
                  end do

                  ! MEANS PRIOR
                  do j = 1,(j0)-1
                        logRho = logRho -0.5d0*(m(j)-xi)**2*kappa
                  end do
                  logRho = logRho &
                              -(1-gt)*0.5d0*(mj0-xi)**2*kappa &
                              -gt*0.5d0*(mj1-xi)**2*kappa &
                              -gt*0.5d0*(mj2-xi)**2*kappa
                  do j = (j0+1)+1,k+1
                        logRho = logRho -0.5d0*(m(j)-xi)**2*kappa
                  end do

                  ! VARIANCES PRIOR
                  do j = 1,(j0)-1
                        logRho = logRho &
                                    -(alpha+1)*log(s(j)) -beta/s(j)
                  end do
                  logRho = logRho &
                              -(1-gt)*(alpha+1)*log(sj0) -(1-gt)*beta/sj0 &
                              -gt*(alpha+1)*log(sj1) -gt*beta/sj1 &
                              -gt*(alpha+1)*log(sj2) -gt*beta/sj2
                  do j = (j0+1)+1,k+1
                        logRho = logRho &
                                    -(alpha+1)*log(s(j)) -beta/s(j)
                  end do

                  ! AUXILIARY VALUES
                  logRho = logRho &
                              +(1-gt)*log(aux1) +(1-gt)*log(1-aux1) &
                              +(1-gt)*log(aux2) +(1-gt)*log(1-aux2)

                  ! JACOBIAN
                  logRho = logRho &
                              -(1-gt)*( &
                                    log(wj0) +abs(mj1-mj2) &
                                    +log(sj1) +log(sj2) &
                                    -log(aux2)-log(1-aux2**2) &
                                    -log(aux3)-log(1-aux3) &
                                    -log(sj0) &
                                    )

                  Snumer = Snumer + exp(logRho)

                  pr(j0) = Snumer

            end do

            ! FINALIZE THE CUMULATIVE PROBABILITIES
            if (Snumer.gt.0.d0) then
                  do j0 = 1,k
                        pr(j0) = pr(j0)/Snumer
                  end do
            else
                  ! OVERFLOW SECURITY
                  do j0 = 1,k
                        pr(j0) = dble(j0)/dble(k)
                  end do
            end if

            ! GENERATE THE LABEL
            call rnguniform(u)
            j0 = 1
            do
                  if (pr(j0).gt.u) exit
                  if (j0.eq.k) exit
                  j0 = j0+1
            end do
            ! CHECK FOR ERRORS
            if (j0.gt.k) stop 'j0 AIS UPDATE::(j0.gt.k)'

      end if

      end do

      end subroutine AisSplitSweep

end module splitmergemove_mod



