      program main

      implicit none
      include 'silo.inc'

      integer ndim, err, dbid, stat, i, j, n, lmname, lvname
      integer dims(3)
      real(4) orig(3), delta(3,3)
      real(4), allocatable :: coords(:,:)
      real(4), allocatable :: wf(:,:,:)
      character(256) fname, bname, mname, vname, pre, ext

      call getarg(1,fname)
      open(80, file=trim(fname), action='READ', form='FORMATTED')
      read(80, '(3I6)')dims(:)
      ndim = max(dims(1), dims(2), dims(3))
      allocate(coords(3,ndim))
      allocate(wf(dims(1),dims(2),dims(3)))
      call readr3d(80, ndim, dims, coords, wf)
      close(80)

      mname = 'quad_mesh'
      lmname = 9
      vname = 'rho'
      lvname = 3

      err = dbcreate('out.silo', 8, 0, DB_LOCAL, 'file info',
     &               9, DB_PDB, dbid)

      err = dbputqm(dbid, mname, lmname, 'x', 1, 'y', 1, 'z', 1,
     &              coords(1,:), coords(2,:), coords(3,:), dims, 3,
     &              DB_FLOAT, DB_COLLINEAR, DB_F77NULL, stat)

      err = dbputqv1(dbid, vname, lvname, mname, lmname, wf,
     &               dims, 3, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT,
     &               DB_F77NULL, stat)

      err = dbclose(dbid)
      deallocate(wf,coords)

      end program
