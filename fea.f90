module fea

    !! This module contains procedures that are common to FEM-analysis

    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial

        !! This subroutine is mainly used to allocate vectors and matrices

        use fedata
        use link1
        use plane42rect

        integer :: maxn, minn, bww, i, e, nen
        integer :: bwtemp=0
        integer, dimension(mdim) :: edof



! Hint for continuum elements:
!        integer, parameter :: mdim = 8
!        integer, dimension(mdim) :: edof

        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn

        ! calculate bw

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
               edof(2*i-1) = 2 * element(e)%ix(i) - 1
               edof(2*i)   = 2 * element(e)%ix(i)
             end do
            maxn=maxval(edof)
            minn=minval(edof)
            bww=maxn-minn+1
            if (bww>=bwtemp) then
                bwtemp=bww
            else
                bwtemp=bwtemp
            end if
            bw=bwtemp
            print *, 'bw='
            print *, bw
        end do





        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            allocate(kmat(bw, neqn))
            print *, 'ERROR in fea/initial'
            print *, 'Band form not implemented -- you need to add your own code here'
            !stop
        end if
        allocate (p(neqn), d(neqn))
        allocate (strain(ne, 3), stress(ne, 3))

        ! Initial stress and strain
        strain = 0
        stress = 0
    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: e
        real(wp), dimension(:), allocatable :: plotval

        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce

        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, p)
        else
            print *, 'ERROR in fea/displ'
            print *, 'Band form not implemented -- you need to add your own code here'
            stop
        end if

        ! Transfer results
        d(1:neqn) = p(1:neqn)

        ! Recover stress
        call recover

        ! Output results
        call output

        ! Plot deformed structure
        call plotmatlabdef('Deformed')

        ! Plot element values
        allocate (plotval(ne))
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                plotval(e) = 0
                print *, 'WARNING in fea/displ: Plot value not set -- you need to add your own code here'
            end if
        end do
        call plotmatlabeval('Stresses',plotval)

    end subroutine displ
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42rect

        integer :: i
        integer :: k
        integer :: j
! Hint for continuum elements:
!        integer, dimension(mdim) :: edof
!        real(wp), dimension(mdim) :: xe
!        real(wp), dimension(mdim) :: re

        ! Build load vector
        p(1:neqn) = 0
        do i = 1, np
            select case(int(loads(i, 1)))
            case( 1 )
            	k=2*(loads(i,2)-1)+loads(i,3)! Build nodal load contribution
            	!print *,k
                p(k) = loads(i,4)
                print*,'p=='
                print*,p
                !print *, 'WARNING in fea/buildload: You need to replace hardcoded nodal load with your code'
            case( 2 )
                !print *,loads
                print *, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                j=loads(i,2)*loads(i,3)
                !print*,j
                p(j)=loads(1,4)


                print*,'p='
                print*,p

            	! Builpd uniformly distributed surface (pressure) load contribution
                !print *, 'ERROR in fea/buildload'
                !print *, 'Distributed loads nooooot defined -- you need to add your own code here'
                !stop
            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do
    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42rect

        integer :: e, i, j,k
        integer :: nen
! Hint for system matrix in band form:
        integer :: irow, icol,band_row
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke
        real(wp), dimension(bw,neqn) :: kband
        !real(wp), dimension(24,24) :: kmat1
! Hint for modal analysis:
!        real(wp), dimension(mdim, mdim) :: me
        real(wp) :: young, area
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, thk
   ! real(wp) :: nu, dens, thk
        ! Reset stiffness matrix
        if (.not. banded) then
            kmat = 0
            ! print *, kmat
        else
            kmat = 0  !!!!!!!!!!

            print*,'ERROR in fea/buildstiff'
            print*,'Band form not implemented -- you need to add your own code here'
            !stop
        end if

        do e = 1, ne

            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
                 xe(2*i-1) = x(element(e)%ix(i),1)
                 xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
                 young = mprop(element(e)%mat)%young
                 nu  = mprop(element(e)%mat)%nu
                 thk  = mprop(element(e)%mat)%thk
                 call plane42rect_ke(xe, young, nu, thk, ke)

                 !print *, 'ERROR in fea/buildstiff:'
                 !print *, 'Stiffness matrix for plane42rect elements not implemented -- you need to add your own code here'
                 ! stop
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do

! Hint: Can you eliminate the loops above by using a different Fortran array syntax?
            else
            do k=1,ne
                  do i = 1,2*nen
                    irow=edof(i)
                    do j=1, 2*nen
                        icol=edof(j)
                        band_row=irow - icol + bw / 2 + 1
                        if (band_row >= 1 .and. band_row <= bw) then  !Check bounds
                        kband(band_row, icol) = kband(band_row, icol) + ke(i, j)
                        end if
                    end do
                  end do
            end do
print*,'##########################'
print *, kband
print*, ne

                print *, 'ERROR in fea/buildstiff'
                print *, 'Band form not implemented -- you need to add our own code here'
                !stop
            end if
        end do

    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, idof,j,band_row
        real(wp) :: penal
        real(wp), dimension(bw,neqn) :: kband
        ! Correct for supports
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)
                end do
            end if
        else
            kband=kmat
            do i = 1, nb
                idof = int(2*(bound(i,1)-1) + bound(i,2))

                do j = max(1, idof - bw / 2), min(neqn, idof + bw / 2)
                    band_row = idof - j + bw / 2 + 1
                    if (band_row >= 1 .and. band_row <= bw) then
                    kband(band_row, j) = 0.0_wp
                    end if
                end do
            end do
            ! set diagonal element to 1
                kband(bw / 2 + 1, idof) = 1.0_wp

print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print*,kband
            print *, 'ERROR in fea/enforce'
            print *, 'Band form not implemented -- you need to add your own code here'
            stop
        end if
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain,
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42rect

        integer :: e, i, nen, eface
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area, nu, thk, fe, re(8)
! Hint for continuum elements:
!        real(wp):: nu, dens, thk
        real(wp), dimension(3) :: estrain, estress

        ! Reset force vector
        p = 0

        do e = 1, ne

            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))   !nen等于4
                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain  !因为estrain中就是有三个元素

                print*, 'stress='
                print*, stress

                print*, 'estress='
                print*, estress


            case( 2 )
                !print *, 'WARNING in fea/recover: Stress and strain not calculated for continuum' &
                !    // 'elements -- you need to add your own code here'



                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!这部分是对于面力的？对点力无效？
                young = mprop(element(e)%mat)%young
                thk  = mprop(element(e)%mat)%thk
                nu  = mprop(element(e)%mat)%nu

                do i=1,SIZE(loads, DIM=1)
                if (loads(i,2)==e) then
                        eface=loads(i,3)
                        fe=loads(i,4)
                        !call plane42rect_re(xe, eface, fe, thk, re)
                        !print*, 'i============'
                        !print*, i
                end if
                   !print*, 'iii============'
                   !print*, i
                end do


                call plane42rect_ke(xe, young, nu, thk, ke)


                call plane42rect_re(xe, eface, fe, thk, re)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))-re


                !print*, 're============'
                !print*, nface




                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                call plane42rect_ss(xe, de, young, nu, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain

                !print*, 'stress=='
                !print*, stress


            end select
        end do
        print*, 'stress============'
        print*, stress
    end subroutine recover

end module fea
