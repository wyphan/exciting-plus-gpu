
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfftifc
! !INTERFACE:
subroutine zfftifc(nd,n,sgn,z)
! !INPUT/OUTPUT PARAMETERS:
!   nd   : number of dimensions (in,integer)
!   n    : grid sizes (in,integer(nd))
!   sgn  : FFT direction, -1: forward; 1: backward (in,integer)
!   z    : array to transform (inout,complex(n(1)*n(2)*...*n(nd)))
! !DESCRIPTION:
!   Interface to the double-precision complex fast Fourier transform routine.
!   This is to allow machine-optimised routines to be used without affecting the
!   rest of the code. See routine {\tt nfftifc}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
  USE mod_prec, ONLY: dz
#ifdef _FFTW3_
  USE mod_lapack, ONLY: ZDSCAL
  USE, INTRINSIC :: ISO_C_BINDING ! For C_PTR and C_DOUBLE_COMPLEX
  implicit none
#elif defined(_MKL_)
  use MKL_DFTI ! this module required by MKL
  IMPLICIT NONE
#else
  USE fftpack5
  IMPLICIT NONE
#endif /* _FFTW3_ */

! arguments
integer, intent(in) :: nd
integer, intent(in) :: sgn
integer, intent(in) :: n(nd)
complex(KIND=dz), INTENT(INOUT), TARGET :: z(*)

#ifdef _FFTW3_
!-------------------------------------!
!     interface to FFTW version 3     !
!-------------------------------------!
INCLUDE 'fftw3.f03'
!integer, parameter :: FFTW_ESTIMATE=64
integer i,p
TYPE(C_PTR) :: plan
real(8) t1
INTEGER(KIND=C_INT) :: nrev(nd)
COMPLEX(KIND=dz), DIMENSION(PRODUCT(n)) :: ztmp

! Reverse grid indices
DO i = 1, nd
   nrev(i) = n(nd-i+1)
END DO

! TODO: Move this into its own subroutine
! TODO: Use wisdom plan if available
!$OMP CRITICAL
plan = fftw_plan_dft( nd, nrev, z, z, sgn, FFTW_ESTIMATE )
!$OMP END CRITICAL

CALL fftw_execute_dft( plan, z, z )

! TODO: Move this into its own subroutine
!$OMP CRITICAL
CALL fftw_destroy_plan(plan)
!$OMP END CRITICAL

! Normalize the result (for backward transform)
if (sgn.eq.-1) then
  p=1
  do i=1,nd
    p=p*n(i)
  end do
  t1=1.d0/dble(p)
  call zdscal(p,t1,z,1)
end if

#elif defined( _MKL_ )
!----------------------------------!
!     interface to MKL 8.1/9.1     !
!----------------------------------!
! (with thanks to Torbjorn Bjorkman)
integer dftistatus,i,p
real(8) dftiscale
type(DFTI_DESCRIPTOR), POINTER :: handle
p=1
do i=1,nd
  p=p*n(i)
end do
dftiscale=1.d0/dble(p)
dftistatus=DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,nd,n)
dftistatus=DftiSetValue(handle, DFTI_FORWARD_SCALE,dftiscale)
dftistatus=DftiCommitDescriptor(handle)
if (sgn.eq.-1) then
  dftistatus=DftiComputeForward(handle,z)
else
  dftistatus=DftiComputeBackward(handle,z)
end if
dftistatus=DftiFreeDescriptor(handle)

#elif defined(_Z3DFFT_)
!-------------------------------------!
!     interface to HP MLIB Z3DFFT     !
!-------------------------------------!
integer ier
if (nd.eq.3) then
  call z3dfft(z,n(1),n(2),n(3),n(1),n(2),sgn,ier)
end if

#else
!----------------------------------------!
!     interface to modified FFTPACK5     !
!----------------------------------------!
call cfftnd(nd,n,sgn,z)
#endif /* _OPENACC */

!-------------------------------------!

return
end subroutine
!EOC

