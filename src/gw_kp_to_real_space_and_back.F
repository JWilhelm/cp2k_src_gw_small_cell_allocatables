!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 05.2024
! **************************************************************************************************
MODULE gw_kp_to_real_space_and_back
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_read_unformatted, cp_fm_release, &
        cp_fm_set_all, cp_fm_to_fm, cp_fm_type, cp_fm_write_unformatted
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig_canon
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_dist2d_to_dist
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_p_type, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE gw_utils,                        ONLY: de_init_bs_env, time_to_freq, analyt_conti_and_print, &
get_VBM_CBM_bandgaps
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: data_3_type,&
                                              post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists
   USE rpa_gw_kpoints_util,             ONLY: cp_cfm_power
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_kp_to_real_space_and_back'

   PUBLIC :: fm_trafo_rs_to_ikp, trafo_rs_to_ikp, trafo_ikp_to_rs, fm_add_ikp_to_rs, &
                  add_ikp_to_all_rs

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cfm_ikp ...
!> \param fm_rs ...
!> \param kpoints ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE fm_trafo_rs_to_ikp(cfm_ikp, fm_rs, kpoints, ikp)
      TYPE(cp_cfm_type)                                  :: cfm_ikp
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_trafo_rs_to_ikp'

      COMPLEX(KIND=dp)                                   :: expikR
      INTEGER                                            :: handle, img, nimages, nimages_fm_rs
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)

      nimages = SIZE(kpoints%index_to_cell, 1)
      nimages_fm_rs = SIZE(fm_rs)

      CPASSERT(nimages == nimages_fm_rs)

      cfm_ikp%local_data(:, :) = z_zero
      DO img = 1, nimages

         CALL add_rs_to_ikp(fm_rs(img)%local_data, cfm_ikp%local_data, kpoints%index_to_cell, &
                            kpoints%xkp(1:3, ikp), img)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_trafo_rs_to_ikp

! **************************************************************************************************
!> \brief ...
!> \param array_rs ...
!> \param array_kp ...
!> \param index_to_cell ...
!> \param xkp ...
! **************************************************************************************************
   SUBROUTINE trafo_rs_to_ikp(array_rs, array_kp, index_to_cell, xkp)
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: array_rs
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      INTEGER, DIMENSION(:, :)                           :: index_to_cell
      REAL(KIND=dp)                                      :: xkp(3)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trafo_rs_to_ikp'

      COMPLEX(KIND=dp)                                   :: expikR
      INTEGER                                            :: handle, i_cell, nimages
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)

      nimages = SIZE(index_to_cell, 1)

      CPASSERT(nimages == SIZE(array_rs, 3))

      array_kp(:, :) = 0.0_dp
      DO i_cell = 1, nimages

         CALL add_rs_to_ikp(array_rs(:, :, i_cell), array_kp, index_to_cell, xkp, i_cell)

      END DO

      CALL timestop(handle)

   END SUBROUTINE trafo_rs_to_ikp

   SUBROUTINE add_rs_to_ikp(array_rs, array_kp, index_to_cell, xkp, i_cell)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      REAL(KIND=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER                                            :: cell(3), i_cell
      REAL(KIND=dp)                                      :: xkp(3)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp
      INTEGER, DIMENSION(:, :)                           :: index_to_cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_rs_to_ikp'

      INTEGER                                            :: handle
      COMPLEX(KIND=dp)                                   :: expikR
      REAL(KIND=dp)                                      :: arg, im, re

      CALL timeset(routineN, handle)

      arg = REAL(index_to_cell(i_cell, 1), dp)*xkp(1) + &
            REAL(index_to_cell(i_cell, 2), dp)*xkp(2) + &
            REAL(index_to_cell(i_cell, 3), dp)*xkp(3)

      expikR = z_one*COS(twopi*arg) + gaussi*SIN(twopi*arg)

      array_kp(:, :) = array_kp(:, :) + expikR*array_rs(:, :)

      CALL timestop(handle)

   END SUBROUTINE add_rs_to_ikp

! **************************************************************************************************
!> \brief ...
!> \param array_kp ...
!> \param array_rs ...
!> \param cell ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE trafo_ikp_to_rs(array_kp, array_rs, cell, kpoints)
      COMPLEX(KIND=dp), DIMENSION(:, :, :)               :: array_kp
      REAL(KIND=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER                                            :: cell(3)
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trafo_ikp_to_rs'

      INTEGER                                            :: handle, ikp

      CALL timeset(routineN, handle)

      CPASSERT(kpoints%nkp == SIZE(array_kp, 3))

      array_rs(:, :) = 0.0_dp

      DO ikp = 1, kpoints%nkp

         CALL add_ikp_to_rs(array_kp(:, :, ikp), array_rs, cell, kpoints, ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE trafo_ikp_to_rs

   SUBROUTINE fm_add_ikp_to_rs(cfm_ikp, fm_rs, kpoints, ikp)
      TYPE(cp_cfm_type)                                  :: cfm_ikp
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_add_ikp_to_rs'

      INTEGER                                            :: handle, img, nimages, nimages_fm_rs
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)

      nimages = SIZE(kpoints%index_to_cell, 1)
      nimages_fm_rs = SIZE(fm_rs)

      CPASSERT(nimages == nimages_fm_rs)

      ALLOCATE (index_to_cell(nimages, 3))
      index_to_cell(1:nimages, 1:3) = kpoints%index_to_cell(1:nimages, 1:3)

      DO img = 1, nimages

         IF (ikp == 1) CALL cp_fm_set_all(fm_rs(img), 0.0_dp)

         CALL add_ikp_to_rs(cfm_ikp%local_data(:, :), fm_rs(img)%local_data, &
                            index_to_cell(img, 1:3), kpoints, ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_add_ikp_to_rs

! **************************************************************************************************
!> \brief ...
!> \param array_kp ...
!> \param array_rs ...
!> \param kpoints ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE add_ikp_to_all_rs(array_kp, array_rs, kpoints, ikp, index_to_cell_ext)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: array_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp
      INTEGER, DIMENSION(:,:), POINTER, OPTIONAL            :: index_to_cell_ext

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_ikp_to_all_rs'

      INTEGER                                            :: cell(3), handle, img, nimages
      INTEGER, DIMENSION(:,:), POINTER :: index_to_cell

      CALL timeset(routineN, handle)

      IF(PRESENT(index_to_cell_ext)) THEN
        index_to_cell => index_to_cell_ext
      ELSE
        index_to_cell => kpoints%index_to_cell
      END IF

      nimages = SIZE(index_to_cell, 1)
      CPASSERT(SIZE(array_rs, 3) == nimages)
      DO img = 1, nimages

         cell(1:3) = index_to_cell(img, 1:3)

         CALL add_ikp_to_rs(array_kp, array_rs(:, :, img), cell, kpoints, ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE add_ikp_to_all_rs

! **************************************************************************************************
!> \brief ...
!> \param array_kp ...
!> \param array_rs ...
!> \param cell ...
!> \param kpoints ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE add_ikp_to_rs(array_kp, array_rs, cell, kpoints, ikp)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: array_kp
      REAL(KIND=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER                                            :: cell(3)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_ikp_to_rs'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: arg, im, re

      CALL timeset(routineN, handle)

      arg = REAL(cell(1), dp)*kpoints%xkp(1, ikp) + &
            REAL(cell(2), dp)*kpoints%xkp(2, ikp) + &
            REAL(cell(3), dp)*kpoints%xkp(3, ikp)

      re = COS(twopi*arg)*kpoints%wkp(ikp)
      im = SIN(twopi*arg)*kpoints%wkp(ikp)

      array_rs(:, :) = array_rs(:, :) + re*REAL(array_kp(:, :)) + im*AIMAG(array_kp(:, :))

      CALL timestop(handle)

   END SUBROUTINE add_ikp_to_rs

END MODULE gw_kp_to_real_space_and_back
