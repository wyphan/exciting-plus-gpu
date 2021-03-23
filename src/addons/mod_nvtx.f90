!-------------------------------------------------------------------------------
! Fortran bindings for a small subset of the NVIDIA Tools Extensions library
! Using the example posted at 
! https://github.com/maxcuda/NVTX_example
MODULE nvtx

  USE ISO_C_BINDING

  PUBLIC :: nvtxEventAttributes
  PUBLIC :: nvtxRangePushA, nvtxRangePushAArgb, nvtxMarkA
  PUBLIC :: nvtxRangePushEx, nvtxRangePop
  PUBLIC :: nvtxStartRange, nvtxEndRange

  integer, private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', Z'00ff0000', Z'00ffffff' ]
  character, private, target :: tempName(256)

  TYPE, BIND(C) :: nvtxEventAttributes
     integer(C_INT16_T) :: version=1
     integer(C_INT16_T) :: size=48 !
     integer(C_INT) :: category=0
     integer(C_INT) :: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT) :: color
     integer(C_INT) :: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT) :: reserved0
     integer(C_INT64_T) :: payload   ! union uint,int,double
     integer(C_INT) :: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII = 1 
     type(C_PTR) :: message  ! ascii char
  END TYPE nvtxEventAttributes

!-------------------------------------------------------------------------------
! Push range
#ifdef _USE_NVTX_
  INTERFACE nvtxRangePush

     ! push range with custom label and standard color
     subroutine nvtxRangePushA(name) bind(C, name="nvtxRangePushA")
       use iso_c_binding, only : c_char
       character(kind=C_CHAR) :: name(256)
     end subroutine nvtxRangePushA

     ! push range with custom label and custom color
     SUBROUTINE nvtxRangePushEx( event ) BIND(C, name="nvtxRangePushEx")
       USE ISO_C_BINDING
       IMPORT :: nvtxEventAttributes
       IMPLICIT NONE
       TYPE(nvtxEventAttributes) :: event
     END SUBROUTINE nvtxRangePushEx

  END INTERFACE nvtxRangePush
#endif /* _USE_NVTX_ */

!-------------------------------------------------------------------------------
! Pop range
#ifdef _USE_NVTX_
  INTERFACE nvtxRangePop
     SUBROUTINE nvtxRangePop() BIND(C, name="nvtxRangePop")
     END SUBROUTINE nvtxRangePop
  END INTERFACE nvtxRangePop
#endif /* _USE_NVTX_ */

!-------------------------------------------------------------------------------
! Place a mark on the timeline with a message
#ifdef _USE_NVTX_
  INTERFACE nvtxMarkA
     subroutine nvtxMarkA(string) bind(C, name="nvtxMarkA")
       use iso_c_binding, only : c_char
       character(kind=c_char) :: string(*)
     end subroutine nvtxMarkA
  END INTERFACE nvtxMarkA
#endif /* _USE_NVTX_ */

!-------------------------------------------------------------------------------
! Name an OS thread
#ifdef _USE_NVTX_
  INTERFACE nvtxNameOsThread
     subroutine nvtxNameOsThread(tid, string) bind(C, name="nvtxNameOsThread")
       use iso_c_binding, only : c_int, c_char
       integer(kind=c_int) :: tid
       character(kind=c_char) :: string(*)
     end subroutine nvtxNameOsThread
  END INTERFACE nvtxNameOsThread
#endif /* _USE_NVTX_ */

!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
! Start range with a predefined color id
  subroutine nvtxStartRange( name, id )
    IMPLICIT NONE

    character(kind=c_char,len=*) :: name
    integer, optional :: id

#ifdef _USE_NVTX_
    type(nvtxEventAttributes) :: event
    CHARACTER(KIND=C_CHAR, LEN=256) :: trimmed_name
    INTEGER :: i

    trimmed_name=trim(name)//c_null_char

  ! move scalar trimmed_name into character array tempName
    do i=1,LEN(trim(name)) + 1
       tempName(i) = trimmed_name(i:i)
    enddo

    if ( .not. present(id) ) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
#endif /* _USE_NVTX_ */

    RETURN
  end subroutine nvtxStartRange

!-------------------------------------------------------------------------------
! End range
  subroutine nvtxEndRange
    IMPLICIT NONE

#ifdef _USE_NVTX_
    call nvtxRangePop
#endif /* _USE_NVTX_ */

    RETURN
  end subroutine nvtxEndRange

!-------------------------------------------------------------------------------

END MODULE nvtx


