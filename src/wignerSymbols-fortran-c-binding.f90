!******************************************************-/
! This source code is subject to the terms of the GNU  -/
! Lesser Public License. If a copy of the LGPL was not -/
! distributed with this file, you can obtain one at    -/
! https://www.gnu.org/licenses/lgpl.html.              -/
! ******************************************************/ 

! ----------------------------------------------------------------
! - Author:       Joey Dumont <joey.dumont@gmail.com>            -
! - Date created: 2013-08-15                                     -
! - Date modded:  2013-08-15                                     -
! - Desc.:        ISO C Binding wrapper for Wigner-3j and        -
! -               -6j symbols. This will allow the use of        -
! -               Fortran subroutines drc3jj, drc3jm and         -
! -               drc6jj.                                        -
! ----------------------------------------------------------------

subroutine drc3jj_wrap(l2, l3, m2, m3, l1min, l1max, thrcof, ndim, ier) bind(C)

  use iso_c_binding
  implicit none

  real(c_double), value, intent(in)           :: l2, l3, m2, m3
  real(c_double), intent(out)                 :: l1min, l1max
  real(c_double), dimension(ndim), intent(out):: thrcof
  integer (c_int), value, intent(in)          :: ndim
  integer (c_int), intent(out)                :: ier

  interface
          SUBROUTINE DRC3JJ (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, IER)
              INTEGER NDIM, IER
              DOUBLE PRECISION L2, L3, M2, M3, L1MIN, L1MAX, THRCOF(NDIM)
          end SUBROUTINE DRC3JJ
          end interface

          call DRC3JJ(l2, l3, m2, m3, l1min, l1max, thrcof, ndim, ier)

end subroutine drc3jj_wrap

subroutine drc6j_wrap(l2, l3, l4, l5, l6, l1min, l1max, sixcof, ndim, ier) bind(C)
      use iso_c_binding
      implicit none

      real(c_double), value, intent(in)           :: l2, l3, l4, l5, l6
      real(c_double), intent(out)                 :: l1min, l1max
      real(c_double), dimension(ndim), intent(out):: sixcof
      integer(c_int), value, intent(in)           :: ndim
      integer(c_int), intent(out)                 :: ier

      interface
          SUBROUTINE DRC6J(L2, L3, L4, L5, L6, L1MIN, L1MAX, SIXCOF, NDIM, IER)
              INTEGER NDIM, IER
              DOUBLE PRECISION L2, L3, L4, L5, L6, L1MIN, L1MAX, SIXCOF(NDIM)
          END SUBROUTINE DRC6J
      end interface

      call DRC6J(l2, l3, l4, l5, l6, l1min, l1max, sixcof, ndim, ier)
end subroutine drc6j_wrap
