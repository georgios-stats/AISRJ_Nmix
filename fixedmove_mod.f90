! 11/09/2010

! INFORMATIONS:
! ======================================================================
! FIXED MOVE 
! LABEL SWICHING        : CI
! BLOCKWISE SCHEME      : SYSTEMATIC SCAN
! ======================================================================
! PRIORS                : DU,DIRICHLET, NORMAL, GAMMA
! ======================================================================

module fixedmove_mod

      private
      public      :: FixedMoveSweep

contains

! MAIN SUBROUTINE: FIXED MOVE ==========================================

      subroutine FixedMoveSweep(k,w,m,s,beta,&
                                    delta,xi,kappa,alpha,ge,ha,&
                                    en,y,&
                                    accPrW,accPrM,accPrS,&
                                    sclw,sclm,scls)

            implicit none

      ! DECLARE THE ARGUMENTS
            integer, intent(in)                 :: k
            integer, intent(in)                 :: en
            double precision, intent(inout)     :: w(k),m(k),s(k)
            double precision, intent(inout)     :: beta
            double precision, intent(in)        :: xi,kappa,alpha,ge,ha,delta
            double precision, intent(in)        :: y(en)
            double precision, intent(out)       :: accPrW,accPrM,accPrS
            double precision, intent(in)        :: sclw,sclm,scls

      ! DECLARE THE LOCAL VARIABLES
            integer           :: j,i
            double precision  :: AA,BB,vec(k),u
            double precision  :: logPostN,logPost,logPropN,logProp
            logical           :: isOrdered

      ! WEIGHT BLOCK ; FORWARD
      ! ======================

      ! COMPUTE THE PROPOSALS
            AA = w(k)
            do j = 1, k-1
                  call rngnormal(vec(j:j),1)
                  vec(j) = w(j)*exp(sclw*vec(j))
                  AA = AA+vec(j)
            end do
            BB = 0.d0
            do j = 1,k-1
                  vec(j) = vec(j)/AA
                  BB = BB +vec(j)
            end do
            vec(k) = 1-BB

      ! COMPUTE THE ACCEPTANCE PROBABILITY

            ! POSTERIOR
            logPost = 0.d0
            logPostN = 0.d0

            ! LIKELIHOOD
            do i = 1,en
                  AA = 0.d0
                  BB = 0.d0
                  do j = 1,k
                        AA = AA +w(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(y(i)-m(j))**2/s(j))
                        BB = BB +vec(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(y(i)-m(j))**2/s(j))
                  end do
                  logPost = logPost+log(AA)
                  logPostN = logPostN+log(BB)
            end do

            ! PRIOR
            do j = 1,k
                  logPost = logPost +(delta-1.d0)*log(w(j))
                  logPostN = logPostN +(delta-1.d0)*log(vec(j))
            end do

            ! PROPOSALS
            logProp = 0.d0
            logPropN = 0.d0
            do j = 1,k
                  logProp = logProp -log(w(j))
                  logPropN = logPropN -log(vec(j))
            end do


            ! ACCEPTANCE PROBABILITY
            accPrW = min(1.d0, &
                        exp(logPostN-logPost+logProp-logPropN) &
                        )

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (accPrW.gt.u) then
                  do j = 1,k
                        w(j) = vec(j)
                  end do
            end if

      ! MEAN BLOCK ; FORWARD
      ! ====================

            ! GENERATE PROPOSALS
            do j = 1,k
                  call rngnormal(vec(j:j),1)
                  vec(j) = m(j)+vec(j)*sclm
            end do

            ! CHECK THE ORDER OF THE PROPOSALS
            isOrdered = .true.
            do j = 1,k-1
                  if (vec(j).gt.vec(j+1)) then
                        isOrdered = .false.
                        exit
                  end if
            end do

            ! COMPUTE THE ACCEPTANCE PROBABILITY

            if (isOrdered) then

                  logPost = 0.d0
                  logPostN = 0.d0

                  ! LIKELIHOOD
                  do i = 1,en
                        AA = 0.d0
                        BB = 0.d0
                        do j = 1,k
                              AA = AA +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(y(i)-m(j))**2/s(j))
                              BB = BB +w(j)/sqrt(s(j)) &
                                          *exp(-0.5d0*(y(i)-vec(j))**2/s(j))
                        end do
                        logPost = logPost+log(AA)
                        logPostN = logPostN+log(BB)
                  end do

                  ! PRIOR
                  do j = 1,k
                        logPost = logPost -0.5d0*(m(j)-xi)**2*kappa
                        logPostN = logPostN -0.5d0*(vec(j)-xi)**2*kappa
                  end do

                  ! PROPOSALS
                  logProp = 0.d0
                  logPropN = 0.d0

                  ! ACCEPTANCE PROBABILITY
                  accPrM = min(1.d0, &
                              exp(logPostN-logPost+logProp-logPropN) &
                              )

            else
                  accPrM = 0.d0
            end if

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (accPrM.gt.u) then
                  do j = 1,k
                        m(j) = vec(j)
                  end do
            end if

      ! VARIANCES BLOCK ; FORWARD
      ! =========================

            ! GENERATE PROPOSALS
            do j = 1,k
                  call rngnormal(vec(j:j),1)
                  vec(j) = s(j)*exp(vec(j)*scls)
            end do

            ! COMPUTE THE ACCEPTANCE PROBABILITY

            logPost = 0.d0
            logPostN = 0.d0
            ! LIKELIHOOD
            do i = 1,en
                  AA = 0.d0
                  BB = 0.d0
                  do j = 1,k
                        AA = AA +w(j)/sqrt(s(j)) &
                                    *exp(-0.5d0*(y(i)-m(j))**2/s(j))
                        BB = BB +w(j)/sqrt(vec(j)) &
                                    *exp(-0.5d0*(y(i)-m(j))**2/vec(j))
                  end do
                  logPost = logPost+log(AA)
                  logPostN = logPostN+log(BB)
            end do

            ! PRIOR
            do j = 1,k
                  logPost = logPost -(alpha+1)*log(s(j))-beta/s(j)
                  logPostN = logPostN -(alpha+1)*log(vec(j))-beta/vec(j)
            end do

            ! PROPOSALS
            logProp = 0.d0
            logPropN = 0.d0
            do j = 1,k
                  logProp = logProp -log(s(j))
                  logPropN = logPropN -log(vec(j))
            end do

            ! ACCEPTANCE PROBABILITY
            accPrS = min(1.d0, &
                        exp(logPostN-logPost+logProp-logPropN) &
                        )

            ! ACCEPT/REJECT
            call rnguniform(u)
            if (accPrS.gt.u) then
                  do j = 1,k
                        s(j) = vec(j)
                  end do
            end if

      ! BETA BLOCK ; FORWARD
      ! ====================

            AA = ge +k*alpha
            BB = ha
            do j = 1,k
                  BB = BB + 1.0d0/s(j)
            end do
            call rnggamma(beta,AA)
            beta = beta/BB

      end subroutine FixedMoveSweep

end module fixedmove_mod




