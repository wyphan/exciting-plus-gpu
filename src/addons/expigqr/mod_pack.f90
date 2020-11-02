MODULE mod_pack

  IMPLICIT NONE

  ! Number of matrix elements in gntuju, set to either 96 or 128 in getmaxgnt()
  INTEGER :: ngntujumax

  ! Maximum angular momentum for local orbitals, set to 3 in gengntuju()
  INTEGER :: lmaxapwlo

  ! total angular momentum lmaxapw + lmaxapwlo
  INTEGER :: Nlmo

END MODULE mod_pack
