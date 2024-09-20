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
MODULE gw_small_cell_full_kp
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE rpa_gw_im_time_util,             ONLY: compute_weight_re_im,&
                                              get_atom_index_from_basis_function_index
 USE gw_kp_to_real_space_and_back, ONLY: fm_trafo_rs_to_ikp, trafo_rs_to_ikp, fm_add_ikp_to_rs, &
                  add_ikp_to_all_rs
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
   USE gw_integrals,                    ONLY: build_3c_integral_block
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_small_cell_full_kp'

   PUBLIC :: gw_calc_small_cell_full_kp

CONTAINS

! **************************************************************************************************
!> \brief Perform GW band structure calculation
!> \param qs_env ...
!> \param bs_env ...
!> \param post_scf_bandstructure_section ...
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE gw_calc_small_cell_full_kp(qs_env, bs_env, post_scf_bandstructure_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_small_cell_full_kp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! G^occ_µλ(i|τ|,k) = sum_n^occ C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
      ! G^vir_µλ(i|τ|,k) = sum_n^vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
      ! k-point k -> cell S: G^occ/vir_µλ^S(i|τ|) = sum_k w_k G^occ/vir_µλ(i|τ|,k) e^(ikS)
      ! χ_PQ^R(iτ) = sum_λR1νR2 [ sum_µS (µR1-S νR2 | P0) G^vir_µλ^S(i|τ|) ]
      !                         [ sum_σS (σR2-S λR1 | QR) G^occ_σν^S(i|τ|) ]
      CALL compute_chi(bs_env, qs_env)

      ! χ_PQ^R(iτ) -> χ_PQ(iω,k) -> ε_PQ(iω,k) -> W_PQ(iω,k) -> Ŵ(iω,k) = M^-1(k)*W(iω,k)*M^-1(k)
      !            -> Ŵ_PQ^R(iτ)
      CALL compute_W_real_space(bs_env, qs_env)

      ! D_µν(k) = sum_n^occ C^*_µn(k) C_νn(k), V^tr_PQ^R = <phi_P,0|V^tr|phi_Q,R>
      ! V^tr(k) = sum_R e^ikR V^tr^R, M(k) = sum_R e^ikR M^R, M(k) -> M^-1(k)
      ! -> Ṽ^tr(k) = M^-1(k) * V^tr(k) * M^-1(k) -> Ṽ^tr_PQ^R = sum_k w_k e^-ikR Ṽ^tr_PQ(k)
      ! Σ^x_λσ^R = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) D_µν^S2    ]
      !                       [ sum_QR2 (σR νS1    | QR1-R2) Ṽ^tr_PQ^R2 ]
      CALL compute_Sigma_x(bs_env, qs_env)

      ! Σ^c_λσ^R(iτ) = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) G^occ/vir_µν^S2(i|τ|) ]
      !                           [ sum_QR2 (σR νS1    | QR1-R2) Ŵ_PQ^R2(iτ)           ]
      CALL compute_Sigma_c(bs_env, qs_env)

      ! Σ^c_λσ^R(iτ,k=0) -> Σ^c_nn(ϵ,k); ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ,k) + Σ^x_nn(k) - v^xc_nn(k)
      CALL compute_QP_energies(bs_env, qs_env)

      CALL de_init_bs_env(bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_small_cell_full_kp

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_chi(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_chi'

      INTEGER :: atom_1, atom_2, atom_P, atom_Q, handle, i_cell_Delta_R, i_cell_R, i_cell_R1, &
         i_cell_R2, i_cell_R_inv, i_t, i_task_local, img, ispin, nimages_3c, P_RI, Q_RI
      INTEGER(KIND=int_8)                                :: n_double_in_3c
      INTEGER, DIMENSION(3)                              :: cell_R, cell_R2, cell_R_inv
      LOGICAL                                            :: is_R1_in_3c_cells, memory_fully_used
      REAL(KIND=dp)                                      :: rP(3), rPQ(3), rQ(3), t1, tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: chi_R, G_occ_S, G_vir_S
      TYPE(cell_type), POINTER                           :: cell
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_occ, M_vir
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      n_double_in_3c = INT(0, KIND=int_8)
      nimages_3c = bs_env%nimages_3c
      ALLOCATE (M_occ(nimages_3c, nimages_3c), M_vir(nimages_3c, nimages_3c))
      ALLOCATE (store_int_3c(nimages_3c, nimages_3c, bs_env%n_atom))
      ! G_µλ^S, S: cell; χ_PQ^R, R: cell, replicated on each MPI process
      ALLOCATE (chi_R(bs_env%n_RI, bs_env%n_RI, bs_env%nimages_scf_desymm))
      ALLOCATE (G_occ_S(bs_env%n_ao, bs_env%n_ao, bs_env%nimages_scf_desymm))
      ALLOCATE (G_vir_S(bs_env%n_ao, bs_env%n_ao, bs_env%nimages_scf_desymm))

      memory_fully_used = .FALSE.

      DO i_t = 1, bs_env%num_time_freq_points

         t1 = m_walltime()
         tau = bs_env%imag_time_points(i_t)

         chi_R(:, :, :) = 0.0_dp

         DO ispin = 1, bs_env%n_spin

            ! 1. compute G^occ,S(iτ) and G^vir^S(iτ) in imaginary time for cell S
            !    Background: G^σ,S(iτ) = G^occ,S,σ(iτ) * Θ(-τ) + G^vir,S,σ(iτ) * Θ(τ), σ ∈ {↑,↓}
            !    G^occ_µλ(i|τ|,k) = sum_n^occ C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
            !    G^vir_µλ(i|τ|,k) = sum_n^vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
            !    k-point k -> cell S: G^occ/vir_µλ^S(i|τ|) = sum_k w_k G^occ/vir_µλ(i|τ|,k) e^(ikS)
            CALL G_occ_vir(bs_env, tau, G_occ_S, ispin, occ=.TRUE., vir=.FALSE.)
            CALL G_occ_vir(bs_env, tau, G_vir_S, ispin, occ=.FALSE., vir=.TRUE.)

            DO i_task_local = 1, bs_env%n_tasks_local

               PRINT *,  'mep =', bs_env%para_env%mepos, &
                  'i_task_local =', i_task_local, ' / ', bs_env%n_tasks_local, &
                  'memused =', memory_fully_used, 'n_d_3c =', n_double_in_3c

               atom_1 = bs_env%task_atom_1_atom_2_img(i_task_local, 1)
               atom_2 = bs_env%task_atom_1_atom_2_img(i_task_local, 2)
               i_cell_Delta_R = bs_env%task_atom_1_atom_2_img(i_task_local, 3)

               CALL reallocate_data_3(M_occ, atom_1, atom_2, i_cell_Delta_R, bs_env, occ=.TRUE.)
               CALL reallocate_data_3(M_vir, atom_1, atom_2, i_cell_Delta_R, bs_env, vir=.TRUE.)

               DO i_cell_R2 = 1, nimages_3c

                  CALL get_second_R(i_cell_Delta_R, i_cell_R2, bs_env, is_R1_in_3c_cells, i_cell_R1)

                  IF (.NOT. is_R1_in_3c_cells) CYCLE

                  ! 2. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator) and
                  !    M_λR1,νR2,P0 = sum_µS (λR1 µR2-S | P0) G_νµ^S for i_task_local
                  CALL compute_3c_and_contract(M_occ, G_occ_S, store_int_3c, atom_1, &
                                               atom_2, i_cell_R1, i_cell_R2, bs_env, qs_env, &
                                               n_double_in_3c, memory_fully_used)
                  CALL compute_3c_and_contract(M_vir, G_vir_S, store_int_3c, atom_2, &
                                               atom_1, i_cell_R2, i_cell_R1, bs_env, qs_env, &
                                               n_double_in_3c, memory_fully_used)

               END DO ! i_cell_R2

               !  3. χ_PQ^R(iτ) = sum_λR1,νR2 M^occ_λR1,νR2,P0 M^vir_νR2,λR1,QR
               CALL contract_M_occ_vir_to_chi(M_occ, M_vir, chi_R, bs_env, i_cell_Delta_R)

            END DO ! i_task_local i.e. i_cell_R1, atom_1, atom_2

         END DO ! ispin

         PRINT *, 'mepos =', bs_env%para_env%mepos, 'χ(iτ =', i_t, &
                  ') exec. time: ', m_walltime() - t1, ' s'

         CALL bs_env%para_env%sync()

         CALL bs_env%para_env%sum(chi_R)

         chi_R(:, :, :) = chi_R(:, :, :)*bs_env%spin_degeneracy

         CALL local_array_to_fm(chi_R, bs_env%fm_chi_R_t(:, i_t))

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',        Execution time', m_walltime() - t1, ' s'
            IF (memory_fully_used .AND. i_t == 1) THEN
               WRITE (bs_env%unit_nr, '(T2,3A)') 'The memory for storing 3-center integrals', &
                  ' is fully used. Increase the computational efficiency by more', &
                  ' memory per core.'
            END IF
         END IF

! JW 2del THIS IS ONLY FOR DEBUGGING!!!!!!!!!!!!!!
         DO i_cell_R = 1, bs_env%nimages_scf_desymm

            cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)

            CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

            DO i_cell_R2 = 1, bs_env%nimages_scf_desymm
               cell_R2 = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R2, 1:3)
               IF (cell_R2(1) == -cell_R(1) .AND. cell_R2(2) == -cell_R(2)) THEN
                  i_cell_R_inv = i_cell_R2
                  cell_R_inv = cell_R2
               END IF
            END DO

            DO P_RI = 1, MIN(bs_env%n_RI, 20)
               DO Q_RI = 1, MIN(bs_env%n_RI, 20)
                  IF (i_t == 1 .AND. bs_env%para_env%mepos == 0) THEN

                     atom_P = (P_RI - 1)/(bs_env%n_RI/bs_env%n_atom) + 1
                     atom_Q = (Q_RI - 1)/(bs_env%n_RI/bs_env%n_atom) + 1
                     rP(1:3) = pbc(particle_set(atom_P)%r(1:3), cell)
                     rQ(1:3) = pbc(particle_set(atom_Q)%r(1:3), cell) + &
                               MATMUL(bs_env%hmat, REAL(cell_R, dp))
                     rPQ(1:3) = rP(1:3) - rQ(1:3)

                     PRINT *, 'P =', P_RI, 'Q =', Q_RI, 'R =', cell_R(1:2), &
                        'rPQ =', rPQ(1:2)*0.529_dp, &
                        'chi_PQ^R =', chi_R(P_RI, Q_RI, i_cell_R)
                     IF (ABS(chi_R(P_RI, Q_RI, i_cell_R)) > 1E-9_dp .AND. &
                         ABS(chi_R(P_RI, Q_RI, i_cell_R) - &
                             chi_R(Q_RI, P_RI, i_cell_R_inv))/ &
                         ABS(chi_R(P_RI, Q_RI, i_cell_R)) > 1E-2_dp) THEN
                        PRINT *, 'ERROR'
                        PRINT *, 'P =', Q_RI, 'Q =', P_RI, 'R =', cell_R_inv(1:2), &
                           'rPQ =', rPQ(1:2)*0.529_dp, &
                           'chi_PQ^R =', chi_R(Q_RI, P_RI, i_cell_R_inv)
                     END IF
                  END IF
               END DO
            END DO

!        DO P_RI = 1, bs_env%n_RI
!           DO Q_RI = 1, bs_env%n_RI
!              el = chi_R(P_RI, Q_RI, i_cell_R)
!              IF(i_t == 1 .AND. bs_env%para_env%mepos == 0 .AND. &
!                 (cell_R(1) .NE. 0 .OR. cell_R(2) .NE. 0) .AND. ABS(el) > 1.0E-3_dp)  THEN
!                 atom_P = (P_RI-1)/(bs_env%n_RI/bs_env%n_atom)+1
!                 atom_Q = (Q_RI-1)/(bs_env%n_RI/bs_env%n_atom)+1
!                 rP(1:3) = pbc(particle_set(atom_P)%r(1:3), cell)
!                 rQ(1:3) = pbc(particle_set(atom_Q)%r(1:3), cell) + &
!                           MATMUL(bs_env%hmat, REAL( cell_R, dp))
!                 rPQ(1:3) = rP(1:3) - rQ(1:3)
!                 PRINT *, 'P =', P_RI, 'Q =', Q_RI, 'R =', cell_R(1:2), &
!                         'rPQ =', rPQ(1:2)*0.529_dp, 'nz chi_PQ^R =', el
!              END IF
!           END DO
!        END DO

         END DO

! END JW 2del

      END DO ! i_t

      CALL timestop(handle)

   END SUBROUTINE compute_chi

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param tau ...
!> \param G_S ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
! **************************************************************************************************
   SUBROUTINE G_occ_vir(bs_env, tau, G_S, ispin, occ, vir)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: G_S
      INTEGER                                            :: ispin
      LOGICAL                                            :: occ, vir

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_occ_vir'

      INTEGER                                            :: handle, homo, ikp, j, j_col_local, n_mo, &
                                                            ncol_local, nimages, nkp
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp)                                      :: tau_E

      CALL timeset(routineN, handle)

      CPASSERT(occ .NEQV. vir)

      CALL cp_cfm_get_info(matrix=bs_env%cfm_work_mo, &
                           ncol_local=ncol_local, &
                           col_indices=col_indices)

      nkp = bs_env%nkp_scf_desymm
      nimages = bs_env%nimages_scf_desymm
      n_mo = bs_env%n_ao
      homo = bs_env%n_occ(ispin)

      DO ikp = 1, nkp

         ! get C_µn(k)
         CALL cp_fm_to_cfm(bs_env%fm_mo_coeff_kp(ikp, ispin, 1), &
                           bs_env%fm_mo_coeff_kp(ikp, ispin, 2), bs_env%cfm_work_mo)

         ! G^occ/vir_µλ(i|τ|,k) = sum_n^occ/vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
         DO j_col_local = 1, ncol_local

            j = col_indices(j_col_local)

            ! 0.5 * |(ϵ_nk-ϵ_F)τ|
            tau_E = ABS(tau*0.5_dp*(bs_env%eigenval_scf(j, ikp, ispin) - bs_env%e_fermi(ispin)))

            IF (tau_E < bs_env%stabilize_exp) THEN
               bs_env%cfm_work_mo%local_data(:, j_col_local) = &
                  bs_env%cfm_work_mo%local_data(:, j_col_local)*EXP(-tau_E)
            ELSE
               bs_env%cfm_work_mo%local_data(:, j_col_local) = z_zero
            END IF

            IF ((occ .AND. j > homo) .OR. (vir .AND. j <= homo)) THEN
               bs_env%cfm_work_mo%local_data(:, j_col_local) = z_zero
            END IF

         END DO

         CALL parallel_gemm(transa="N", transb="C", m=n_mo, n=n_mo, k=n_mo, alpha=z_one, &
                            matrix_a=bs_env%cfm_work_mo, matrix_b=bs_env%cfm_work_mo, &
                            beta=z_zero, matrix_c=bs_env%cfm_work_mo_2)

         ! trafo k-point k -> cell S:  G^occ/vir_µλ(i|τ|,k) -> G^occ/vir,S_µλ(i|τ|)
         CALL fm_add_ikp_to_rs(bs_env%cfm_work_mo_2, bs_env%fm_G_S, &
                                    bs_env%kpoints_scf_desymm, ikp)

      END DO ! ikp

      CALL fm_to_local_array(bs_env%fm_G_S, G_S)

      CALL timestop(handle)

   END SUBROUTINE G_occ_vir

! **************************************************************************************************
!> \brief ...
!> \param fm_S ...
!> \param array_S ...
! **************************************************************************************************
   SUBROUTINE fm_to_local_array(fm_S, array_S, weight, add)

      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_S
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: array_S, array_tmp
      REAL(KIND=dp), OPTIONAL :: weight
      LOGICAL, OPTIONAL                                  :: add

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_to_local_array'

      INTEGER                                            :: handle, i, i_row_local, img, j, &
                                                            j_col_local, n_basis, ncol_local, &
                                                            nimages, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_add
      REAL(KIND=dp)                                      :: my_weight

      CALL timeset(routineN, handle)

      my_weight = 1.0_dp
      IF (PRESENT(weight)) my_weight = weight

      my_add = .FALSE.
      IF (PRESENT(add)) my_add = add

      n_basis = SIZE(array_S, 1)
      nimages = SIZE(array_S, 3)

      ! checks
      CPASSERT(SIZE(array_S, 2) == n_basis)
      CPASSERT(SIZE(fm_S) == nimages)
      CPASSERT(LBOUND(array_S, 1) == 1)
      CPASSERT(LBOUND(array_S, 2) == 1)
      CPASSERT(LBOUND(array_S, 3) == 1)

      CALL cp_fm_get_info(matrix=fm_S(1), &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      IF(.NOT. my_add) array_S(:, :, :) = 0.0_dp
      ALLOCATE(array_tmp(SIZE(array_S, 1), SIZE(array_S, 2), SIZE(array_S, 3)))
      array_tmp(:,:,:) = 0.0_dp

      DO img = 1, nimages
         DO i_row_local = 1, nrow_local

            i = row_indices(i_row_local)

            DO j_col_local = 1, ncol_local

               j = col_indices(j_col_local)

               array_tmp(i, j, img) = fm_S(img)%local_data(i_row_local, j_col_local)

            END DO ! j_col_local
         END DO ! i_row_local
      END DO ! img

      CALL fm_S(1)%matrix_struct%para_env%sync()
      CALL fm_S(1)%matrix_struct%para_env%sum(array_tmp)
      CALL fm_S(1)%matrix_struct%para_env%sync()

      array_S(:, :, :) = array_S(:, :, :) + my_weight*array_tmp(:,:,:)

      CALL timestop(handle)

   END SUBROUTINE fm_to_local_array

! **************************************************************************************************
!> \brief ...
!> \param array_S ...
!> \param fm_S ...
!> \param weight ...
!> \param add ...
! **************************************************************************************************
   SUBROUTINE local_array_to_fm(array_S, fm_S, weight, add)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: array_S
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_S
      REAL(KIND=dp), OPTIONAL                            :: weight
      LOGICAL, OPTIONAL                                  :: add

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'local_array_to_fm'

      INTEGER                                            :: handle, i, i_row_local, img, j, &
                                                            j_col_local, n_basis, ncol_local, &
                                                            nimages, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_add
      REAL(KIND=dp)                                      :: my_weight, S_ij

      CALL timeset(routineN, handle)

      my_weight = 1.0_dp
      IF (PRESENT(weight)) my_weight = weight

      my_add = .FALSE.
      IF (PRESENT(add)) my_add = add

      n_basis = SIZE(array_S, 1)
      nimages = SIZE(array_S, 3)

      ! checks
      CPASSERT(SIZE(array_S, 2) == n_basis)
      CPASSERT(SIZE(fm_S) == nimages)
      CPASSERT(LBOUND(array_S, 1) == 1)
      CPASSERT(LBOUND(array_S, 2) == 1)
      CPASSERT(LBOUND(array_S, 3) == 1)

      CALL cp_fm_get_info(matrix=fm_S(1), &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO img = 1, nimages

         DO i_row_local = 1, nrow_local

            i = row_indices(i_row_local)

            DO j_col_local = 1, ncol_local

               j = col_indices(j_col_local)

               IF (my_add) THEN
                  S_ij = fm_S(img)%local_data(i_row_local, j_col_local) + &
                         array_S(i, j, img)*my_weight
               ELSE
                  S_ij = array_S(i, j, img)*my_weight
               END IF
               fm_S(img)%local_data(i_row_local, j_col_local) = S_ij

            END DO ! j_col_local

         END DO ! i_row_local

      END DO ! img

      CALL timestop(handle)

   END SUBROUTINE local_array_to_fm

! **************************************************************************************************
!> \brief ...
!> \param i_cell_Delta_R ...
!> \param i_cell_R2 ...
!> \param bs_env ...
!> \param is_R1_in_3c_cells ...
!> \param i_cell_R1 ...
! **************************************************************************************************
   SUBROUTINE get_second_R(i_cell_Delta_R, i_cell_R2, bs_env, is_R1_in_3c_cells, i_cell_R1)

      INTEGER                                            :: i_cell_Delta_R, i_cell_R2
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL                                            :: is_R1_in_3c_cells
      INTEGER                                            :: i_cell_R1

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_second_R'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(3)                              :: cell_Delta_R, cell_R1, cell_R2

      CALL timeset(routineN, handle)

      cell_R2(1:3) = bs_env%index_to_cell_3c(i_cell_R2, 1:3)

      cell_Delta_R(1:3) = bs_env%index_to_cell_Delta_R(i_cell_Delta_R, 1:3)

      cell_R1(1:3) = cell_R2(1:3) + cell_Delta_R(1:3)

      CALL is_cell_in_index_to_cell(cell_R1, bs_env%index_to_cell_3c, is_R1_in_3c_cells)

      IF (is_R1_in_3c_cells) THEN
         i_cell_R1 = bs_env%cell_to_index_3c(cell_R1(1), cell_R1(2), cell_R1(3))
      ELSE
         i_cell_R1 = -1000
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_second_R

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param index_to_cell ...
!> \param cell_found ...
! **************************************************************************************************
   SUBROUTINE is_cell_in_index_to_cell(cell, index_to_cell, cell_found)
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :)                           :: index_to_cell
      LOGICAL                                            :: cell_found

      CHARACTER(LEN=*), PARAMETER :: routineN = 'is_cell_in_index_to_cell'

      INTEGER                                            :: handle, i_cell, nimg
      INTEGER, DIMENSION(3)                              :: cell_i

      CALL timeset(routineN, handle)

      nimg = SIZE(index_to_cell, 1)

      cell_found = .FALSE.

      DO i_cell = 1, nimg

         cell_i(1:3) = index_to_cell(i_cell, 1:3)

         IF (cell_i(1) == cell(1) .AND. cell_i(2) == cell(2) .AND. cell_i(3) == cell(3)) THEN
            cell_found = .TRUE.
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE is_cell_in_index_to_cell

! **************************************************************************************************
!> \brief ...
!> \param M ...
!> \param atom_1 ...
!> \param atom_2 ...
!> \param i_cell_Delta_R ...
!> \param bs_env ...
!> \param occ ...
!> \param vir ...
! **************************************************************************************************
   SUBROUTINE reallocate_data_3(M, atom_1, atom_2, i_cell_Delta_R, bs_env, occ, vir)
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M
      INTEGER                                            :: atom_1, atom_2, i_cell_Delta_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL, OPTIONAL                                  :: occ, vir

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reallocate_data_3'

      INTEGER                                            :: handle, i_ao_1_end, i_ao_1_size, &
                                                            i_ao_1_start, i_ao_2_end, i_ao_2_size, &
                                                            i_ao_2_start, i_cell_R1, i_cell_R2
      LOGICAL                                            :: is_R1_in_3c_cells, my_occ, my_vir

      IF (PRESENT(occ)) THEN
         my_occ = occ
         my_vir = .NOT. occ
         CPASSERT(.NOT. PRESENT(vir))
      END IF
      IF (PRESENT(vir)) THEN
         my_vir = vir
         my_occ = .NOT. vir
         CPASSERT(.NOT. PRESENT(occ))
      END IF
      CPASSERT(PRESENT(occ) .OR. PRESENT(vir))

      CALL timeset(routineN, handle)

      i_ao_1_start = bs_env%i_ao_start_from_atom(atom_1)
      i_ao_2_start = bs_env%i_ao_start_from_atom(atom_2)
      i_ao_1_end = bs_env%i_ao_end_from_atom(atom_1)
      i_ao_2_end = bs_env%i_ao_end_from_atom(atom_2)

      i_ao_1_size = i_ao_1_end - i_ao_1_start + 1
      i_ao_2_size = i_ao_2_end - i_ao_2_start + 1

      DO i_cell_R1 = 1, bs_env%nimages_3c
         DO i_cell_R2 = 1, bs_env%nimages_3c

            IF (ALLOCATED(M(i_cell_R1, i_cell_R2)%data_3)) THEN
               DEALLOCATE (M(i_cell_R1, i_cell_R2)%data_3)
            END IF

         END DO
      END DO

      DO i_cell_R2 = 1, bs_env%nimages_3c

         CALL get_second_R(i_cell_Delta_R, i_cell_R2, bs_env, is_R1_in_3c_cells, i_cell_R1)

         IF (.NOT. is_R1_in_3c_cells) CYCLE

         IF (my_occ) THEN
            ALLOCATE (M(i_cell_R1, i_cell_R2)%data_3(1:i_ao_1_size, 1:i_ao_2_size, 1:bs_env%n_RI))
            M(i_cell_R1, i_cell_R2)%data_3(:, :, :) = 0.0_dp
         END IF
         IF (my_vir) THEN
            ALLOCATE (M(i_cell_R2, i_cell_R1)%data_3(1:i_ao_2_size, 1:i_ao_1_size, 1:bs_env%n_RI))
            M(i_cell_R2, i_cell_R1)%data_3(:, :, :) = 0.0_dp
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE reallocate_data_3

! **************************************************************************************************
!> \brief ...
!> \param M ...
!> \param G_S ...
!> \param store_int_3c ...
!> \param atom_1 ...
!> \param atom_2 ...
!> \param i_cell_R1 ...
!> \param i_cell_R2 ...
!> \param bs_env ...
!> \param qs_env ...
!> \param n_double_in_3c ...
!> \param memory_fully_used ...
! **************************************************************************************************
   SUBROUTINE compute_3c_and_contract(M, G_S, store_int_3c, atom_1, atom_2, &
                                      i_cell_R1, i_cell_R2, bs_env, qs_env, &
                                      n_double_in_3c, memory_fully_used)

      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: G_S
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c
      INTEGER                                            :: atom_1, atom_2, i_cell_R1, i_cell_R2
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER(KIND=int_8)                                :: n_double_in_3c
      LOGICAL                                            :: memory_fully_used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_and_contract'

      INTEGER :: handle, handle2, i_ao_end_1, i_ao_end_2, i_ao_size_1, i_ao_size_2, i_ao_start_1, &
         i_ao_start_2, i_cell_j, i_cell_k, i_cell_S, i_mu, i_nu, n_ao
      INTEGER, DIMENSION(3)                              :: cell_R1, cell_R1_plus_cell_S, cell_R2, &
                                                            cell_S
      LOGICAL                                            :: cell_found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c, M_tmp

      CALL timeset(routineN, handle)

      cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
      cell_R2(1:3) = bs_env%index_to_cell_3c(i_cell_R2, 1:3)

      n_ao = bs_env%n_ao

      i_ao_start_1 = bs_env%i_ao_start_from_atom(atom_1)
      i_ao_end_1 = bs_env%i_ao_end_from_atom(atom_1)
      i_ao_size_1 = i_ao_end_1 - i_ao_start_1 + 1

      ALLOCATE (int_3c(1:i_ao_size_1, 1:n_ao, 1:bs_env%n_RI))

      i_ao_start_2 = bs_env%i_ao_start_from_atom(atom_2)
      i_ao_end_2 = bs_env%i_ao_end_from_atom(atom_2)
      i_ao_size_2 = i_ao_end_2 - i_ao_start_2 + 1

      ALLOCATE (M_tmp(1:i_ao_size_1, 1:i_ao_size_2, 1:bs_env%n_RI))
      M_tmp(:, :, :) = 0.0_dp

      DO i_cell_S = 1, bs_env%nimages_scf_desymm

         cell_S(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_S, 1:3)
         cell_R1_plus_cell_S(1:3) = cell_R1(1:3) + cell_S(1:3)

         CALL is_cell_in_index_to_cell(cell_R1_plus_cell_S, bs_env%index_to_cell_3c, cell_found)

         IF (.NOT. cell_found) CYCLE

         i_cell_j = bs_env%cell_to_index_3c(cell_R2(1), cell_R2(2), cell_R2(3))
         i_cell_k = bs_env%cell_to_index_3c(cell_R1_plus_cell_S(1), cell_R1_plus_cell_S(2), &
                                            cell_R1_plus_cell_S(3))

         IF(bs_env%nblocks_3c(i_cell_j, i_cell_k) == 0) CYCLE

         IF (ALLOCATED(store_int_3c(i_cell_j, i_cell_k, atom_1)%data_3)) THEN
            int_3c(:, :, :) = store_int_3c(i_cell_j, i_cell_k, atom_1)%data_3(:, :, :)
         ELSE
            ! compute 3-c integrals ( µ (atom_1, cell_j) ν (all at., cell_k) | P(all at., cell 0) )
            ! ("|": truncated Coulomb operator), inside build_3c_integral_block: (j k | i)
            CALL build_3c_integral_block(int_3c, qs_env, bs_env%ri_metric, &
                                    basis_j=bs_env%basis_set_AO, &
                                    basis_k=bs_env%basis_set_AO, &
                                    basis_i=bs_env%basis_set_RI, &
                                    cell_j=cell_R2, &
                                    cell_k=cell_R1_plus_cell_S, &
                                    atom_j=atom_1, &
                                    k_bf_start_from_atom=bs_env%i_ao_start_from_atom, &
                                    i_bf_start_from_atom=bs_env%i_RI_start_from_atom)

            IF (n_double_in_3c < bs_env%avail_doubles_per_proc) THEN
               ALLOCATE (store_int_3c(i_cell_j, i_cell_k, atom_1)% &
                         data_3(1:i_ao_size_1, 1:n_ao, 1:bs_env%n_RI))
               store_int_3c(i_cell_j, i_cell_k, atom_1)%data_3(:, :, :) = int_3c(:, :, :)
               n_double_in_3c = n_double_in_3c + &
                                INT(i_ao_size_1*n_ao, KIND=int_8)*INT(bs_env%n_RI, KIND=int_8)
            ELSE
               ! print warning when memory is full for first time
               IF(.NOT. memory_fully_used) CALL warning(bs_env)
               memory_fully_used = .TRUE.
            END IF
         END IF

         CALL timeset(routineN//"_3c_x_G", handle2)

         ! M_λR1,νR2,P0 = sum_µS (λR1 µR2-S | P0) G_νµ^S for i_task_local
         DO i_nu = 1, i_ao_size_2
            DO i_mu = 1, n_ao
               M_tmp(:, i_nu, :) = M_tmp(:, i_nu, :) + &
                                   int_3c(:, i_mu, :)*G_S(i_nu + i_ao_start_2 - 1, i_mu, i_cell_S)
            END DO
         END DO

         CALL timestop(handle2)

      END DO ! i_cell_S

      M(i_cell_R1, i_cell_R2)%data_3(:, :, :) = M(i_cell_R1, i_cell_R2)%data_3(:, :, :) &
                                                + M_tmp(:, :, :)

      CALL timestop(handle)

   END SUBROUTINE compute_3c_and_contract

! **************************************************************************************************
!> \brief ...
!> \param M_occ ...
!> \param M_vir ...
!> \param chi_R ...
!> \param bs_env ...
!> \param i_cell_Delta_R ...
! **************************************************************************************************
   SUBROUTINE contract_M_occ_vir_to_chi(M_occ, M_vir, chi_R, bs_env, i_cell_Delta_R)
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_occ, M_vir
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: chi_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: i_cell_Delta_R

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_M_occ_vir_to_chi'

      INTEGER :: handle, i_ao_1, i_ao_2, i_cell_R, i_cell_R1, i_cell_R1_minus_R, i_cell_R2, &
         i_cell_R2_minus_R, n_ao_1, n_ao_2, nimages_3c, P_RI
      INTEGER, DIMENSION(3)                              :: cell_R
      LOGICAL                                            :: is_R1_in_3c_cells, &
                                                            is_R1_minus_R_in_3c_cells, &
                                                            is_R2_minus_R_in_3c_cells
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_RI

      CALL timeset(routineN, handle)

      nimages_3c = bs_env%nimages_3c

      ALLOCATE (vec_RI(bs_env%n_RI))

      ! χ_PQ^R(iτ) = sum_λR1,νR2 M^occ_λR1,νR2,P0 M^vir_νR2,λR1,QR
      DO i_cell_R = 1, bs_env%nimages_scf_desymm

         cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)

         DO i_cell_R2 = 1, nimages_3c

            CALL get_second_R(i_cell_Delta_R, i_cell_R2, bs_env, is_R1_in_3c_cells, i_cell_R1)

            IF (.NOT. is_R1_in_3c_cells) CYCLE

            CALL get_Delta_R(i_cell_R1_minus_R, is_R1_minus_R_in_3c_cells, &
                             i_cell_R1, i_cell_R, bs_env)
            CALL get_Delta_R(i_cell_R2_minus_R, is_R2_minus_R_in_3c_cells, &
                             i_cell_R2, i_cell_R, bs_env)

            IF (.NOT. is_R1_minus_R_in_3c_cells) CYCLE
            IF (.NOT. is_R2_minus_R_in_3c_cells) CYCLE

            n_ao_1 = SIZE(M_occ(i_cell_R1, i_cell_R2)%data_3, 1)
            n_ao_2 = SIZE(M_occ(i_cell_R1, i_cell_R2)%data_3, 2)

            DO P_RI = 1, bs_env%n_RI

               vec_RI(:) = 0.0_dp

               DO i_ao_1 = 1, n_ao_1

                  DO i_ao_2 = 1, n_ao_2

                     vec_RI(:) = vec_RI(:) + &
                                 M_occ(i_cell_R1, i_cell_R2)% &
                                 data_3(i_ao_1, i_ao_2, P_RI)* &
                                 M_vir(i_cell_R2_minus_R, i_cell_R1_minus_R)% &
                                 data_3(i_ao_2, i_ao_1, :)
                  END DO

               END DO

               chi_R(P_RI, :, i_cell_R) = chi_R(P_RI, :, i_cell_R) + vec_RI(:)

            END DO ! P_RI

         END DO ! i_ell_R2

      END DO ! i_cell_R

      CALL timestop(handle)

   END SUBROUTINE contract_M_occ_vir_to_chi

! **************************************************************************************************
!> \brief ...
!> \param i_cell_R1_minus_R ...
!> \param cell_found ...
!> \param i_cell_R1 ...
!> \param i_cell_R ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE get_Delta_R(i_cell_R1_minus_R, cell_found, i_cell_R1, i_cell_R, bs_env)
      INTEGER                                            :: i_cell_R1_minus_R
      LOGICAL                                            :: cell_found
      INTEGER                                            :: i_cell_R1, i_cell_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_Delta_R'

      INTEGER                                            :: handle, i_cell
      INTEGER, DIMENSION(3)                              :: cell_i, cell_R, cell_R1, cell_R1_minus_R

      CALL timeset(routineN, handle)

      cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)

      cell_R(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)

      cell_R1_minus_R(1:3) = cell_R1(1:3) - cell_R(1:3)

      cell_found = .FALSE.

      DO i_cell = 1, bs_env%nimages_3c

         cell_i(1:3) = bs_env%index_to_cell_3c(i_cell, 1:3)

         IF (cell_i(1) == cell_R1_minus_R(1) .AND. cell_i(2) == cell_R1_minus_R(2) .AND. &
             cell_i(3) == cell_R1_minus_R(3)) THEN

            cell_found = .TRUE.

            i_cell_R1_minus_R = i_cell

         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE get_Delta_R

   SUBROUTINE warning(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(len=1000)             :: msg

      WRITE (UNIT=msg, FMT=*) &
             "The available memory "//&
             "is too small for storing all 3-center integrals. CP2K therefore recomputes "//&
             "3-center integrals which slows down the calculation substantially. "//&
             "Possible solutions to "//&
             "increase computational efficiency: 1) Use more nodes. 2) Use nodes with "//&
             "more memory or more OMP threads and less MPI ranks (then you need to "//&
             "increase MEMORY_PER_PROC)."

      CPWARN(TRIM(msg))

   END SUBROUTINE warning

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_W_real_space(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_W_real_space'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chi_k_w, eps_k_w, W_k_w, work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: M_inv, M_inv_V_sqrt, V_sqrt, &
                 V_tr_BvK_sqrt, M_inv_V_tr_BvK_sqrt
      INTEGER                                            :: handle, i_cell, i_cell_R, i_cell_R2, &
                                                            i_cell_R_inv, i_t, ikp, ikp_local, &
                                                            j_w, n_RI, nkp, P_RI, Q_RI, i_P, j_Q, &
                                                            nimages_scf_desymm
      INTEGER, DIMENSION(3)                              :: cell_R, cell_R2, cell_R_inv
      REAL(KIND=dp)                                      :: cell(3), el, freq_j, kR, time_i, &
                                                            weight_ij, wkp, xkp(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: chi_R, W_R, MWM_R

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      nimages_scf_desymm = bs_env%nimages_scf_desymm

      ALLOCATE (chi_k_w(n_RI, n_RI), work(n_RI, n_RI), eps_k_w(n_RI, n_RI), W_k_w(n_RI, n_RI))
      ALLOCATE (chi_R(n_RI, n_RI, nimages_scf_desymm), W_R(n_RI, n_RI, nimages_scf_desymm), &
                MWM_R(n_RI, n_RI, nimages_scf_desymm))

      CALL compute_Minv_and_Vsqrt(bs_env, qs_env, M_inv_V_sqrt, M_inv, V_sqrt, V_tr_BvK_sqrt, &
                                  M_inv_V_tr_BvK_sqrt)

      DO j_w = 1, bs_env%num_time_freq_points

         ! χ_PQ^R(iτ) -> χ_PQ^R(iω_j) (which is stored in chi_R, single ω_j from j_w loop)
         chi_R(:, :, :) = 0.0_dp
         DO i_t = 1, bs_env%num_time_freq_points
            freq_j = bs_env%imag_freq_points(j_w)
            time_i = bs_env%imag_time_points(i_t)
            weight_ij = bs_env%weights_cos_t_to_w(j_w, i_t)*COS(time_i*freq_j)

            CALL fm_to_local_array(bs_env%fm_chi_R_t(:, i_t), chi_R, weight_ij, add=.TRUE.)
         END DO

         ikp_local = 0
         W_R(:, :, :) = 0.0_dp
         DO ikp = 1, bs_env%nkp_chi_eps_W_orig_plus_extra

            ! trivial parallelization over k-points
            IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

            ikp_local = ikp_local + 1

            ! 1. χ_PQ^R(iω_j) -> χ_PQ(iω_j,k)
            CALL trafo_rs_to_ikp(chi_R, chi_k_w, bs_env%kpoints_scf_desymm%index_to_cell, &
                                 bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))

            ! JW 2del
            IF (ikp == 1 .AND. j_w == 1) THEN
                 DO i_P = 1, MIN(20,bs_env%n_RI)
                 DO j_Q = 1, MIN(20,bs_env%n_RI)
                  PRINT *, 'chi_PQ(ω_1) at xkp = ', bs_env%kpoints_chi_eps_W%xkp(1:3,1), &
                         i_P, j_Q, chi_k_w(i_P, j_Q)
                 END DO
                 END DO
            END IF

            ! 2. remove negative eigenvalues from χ_PQ(iω,k)
            CALL power(chi_k_w, 1.0_dp, bs_env%eps_eigval_mat_RI)

            ! 3. ε(iω_j,k_i) = Id - V^0.5(k_i)*M^-1(k_i)*χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)

            ! 3. a) work = χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)
            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, chi_k_w, n_RI, &
                       M_inv_V_sqrt(:, :, ikp_local), n_RI, z_zero, work, n_RI)

            ! 3. b) eps_work = V^0.5(k_i)*M^-1(k_i)*work
            CALL ZGEMM('C', 'N', n_RI, n_RI, n_RI, z_one, M_inv_V_sqrt(:, :, ikp_local), n_RI, &
                       work, n_RI, z_zero, eps_k_w, n_RI)

            ! 3. c) ε(iω_j,k_i) = eps_work - Id
            CALL add_on_diag(eps_k_w, z_one)

            ! 4. W(iω_j,k_i) = M^-1(k_i)*V^0.5(k_i)*(ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)*M^-1(k_i)

            ! 4. a) Inversion of ε(iω_j,k_i) using its Cholesky decomposition
            CALL power(eps_k_w, -1.0_dp, 0.0_dp)

            ! 4. b) ε^-1(iω_j,k_i)-Id
            CALL add_on_diag(eps_k_w, -z_one)

            ! JW 2del
            IF (ikp == 1 .AND. j_w == 1) THEN
                 DO i_P = 1, MIN(20,bs_env%n_RI)
                 DO j_Q = 1, MIN(20,bs_env%n_RI)
                  PRINT *, 'ε_PQ(ω_1) at xkp = ', bs_env%kpoints_chi_eps_W%xkp(1:3,1), &
                         i_P, j_Q, eps_k_w(i_P, j_Q)
                 END DO
                 END DO
            END IF

! version 1 JW HACK: INCLUDE V^0.5(k_i)
            ! 4. c) work = (ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)
            CALL ZGEMM('N', 'C', n_RI, n_RI, n_RI, z_one, eps_k_w, n_RI, &
                       V_sqrt(:,:,ikp_local), n_RI, z_zero, work, n_RI)

            ! 4. d) W(iω,k_i) = V^0.5(k_i)*work
            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, V_sqrt(:,:,ikp_local), n_RI, &
                       work, n_RI, z_zero, W_k_w, n_RI)

            ! JW 2del
            IF (ikp == 1 .AND. j_w == 1) THEN
                 DO i_P = 1, MIN(20,bs_env%n_RI)
                 DO j_Q = 1, MIN(20,bs_env%n_RI)
                  PRINT *, 'W_PQ(ω_1) = V^0.5*ε*V^0.5 at xkp = ', &
                         bs_env%kpoints_chi_eps_W%xkp(1:3,1), &
                         i_P, j_Q, W_k_w(i_P, j_Q)
                 END DO
                 END DO
            END IF

! version 2 JW HACK: INCLUDE V^0.5(k_i)*M^-1(k_i)
!
!            ! 4. c) work = (ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)*M^-1(k_i)
!            CALL ZGEMM('N', 'C', n_RI, n_RI, n_RI, z_one, eps_k_w, n_RI, &
!                       M_inv_V_sqrt(:, :, ikp_local), n_RI, z_zero, work, n_RI)
!
!            ! 4. d) W(iω,k_i) = M^-1(k_i)*V^0.5(k_i)*work
!            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, M_inv_V_sqrt(:, :, ikp_local), n_RI, &
!                       work, n_RI, z_zero, W_k_w, n_RI)
!
! version 3 JW HACK: INCLUDE V^tr,BvK^0.5(k_i)*M^-1(k_i)
!
!            ! 4. c) work = (ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)*M^-1(k_i)
!            CALL ZGEMM('N', 'C', n_RI, n_RI, n_RI, z_one, eps_k_w, n_RI, &
!                       M_inv_V_tr_BvK_sqrt(:, :, ikp_local), n_RI, z_zero, work, n_RI)
!
!            ! 4. d) W(iω,k_i) = M^-1(k_i)*V^0.5(k_i)*work
!            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, M_inv_V_tr_BvK_sqrt(:, :, ikp_local), &
!                       n_RI, work, n_RI, z_zero, W_k_w, n_RI)

! version 4 JW HACK: INCLUDE V^tr,BvK^0.5(k_i)
!
!            ! 4. c) work = (ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)
!            CALL ZGEMM('N', 'C', n_RI, n_RI, n_RI, z_one, eps_k_w, n_RI, &
!                       V_tr_BvK_sqrt(:, :, ikp_local), n_RI, z_zero, work, n_RI)
!
!            ! 4. d) W(iω,k_i) = V^0.5(k_i)*work
!            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, V_tr_BvK_sqrt(:, :, ikp_local), &
!                       n_RI, work, n_RI, z_zero, W_k_w, n_RI)

            ! 5. W(iω,k_i) -> W^R(iω) = sum_k w_k e^(-ikR) W(iω,k) [k-point extrapolation here]
            CALL add_ikp_to_all_rs(W_k_w, W_R, bs_env%kpoints_chi_eps_W, ikp, &
                                   index_to_cell_ext=bs_env%kpoints_scf_desymm%index_to_cell)

         END DO ! ikp

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(W_R)

! version 1 JW HACK: INCLUDE M^-1(k) via 
         ! W(iω,k) -> W^R(iω) -> W^R^trunc(iω) -> W^trunc(iω,k) 
         !         -> Ŵ(iω,k) = M^-1(k)*W^trunc(iω,k)*M^-1(k) -> Ŵ^R(iω) (stored in MWM_R)
         CALL truncate_W_and_mult_W_with_Minv(W_R, MWM_R, bs_env, qs_env, j_w)
! end version 1

         ! 6. W^R(iω) -> W^R(iτ) and to fully distributed fm matrix bs_env%fm_MWM_R_t
         DO i_t = 1, bs_env%num_time_freq_points

            freq_j = bs_env%imag_freq_points(j_w)
            time_i = bs_env%imag_time_points(i_t)
            weight_ij = bs_env%weights_cos_w_to_t(i_t, j_w)*COS(time_i*freq_j)
            CALL local_array_to_fm(MWM_R, bs_env%fm_MWM_R_t(:, i_t), weight_ij, add=.TRUE.)
! JW HACK
            CALL local_array_to_fm(W_R, bs_env%fm_W_R_t(:, i_t), weight_ij, add=.TRUE.)
         END DO ! i_t

      END DO ! j_w

! JW 2del THIS IS ONLY FOR DEBUGGING!!!!!!!!!!!!!!
!      CALL fm_to_local_array(bs_env%fm_MWM_R_t(:,1), W_R)
      DO i_t = 1, bs_env%num_time_freq_points

         IF (bs_env%para_env%mepos == 0) PRINT *, 'i_t =', i_t, &
                   'tau =', bs_env%imag_time_points(i_t)

         CALL fm_to_local_array(bs_env%fm_MWM_R_t(:, i_t), MWM_R)
         CALL fm_to_local_array(bs_env%fm_W_R_t(:, i_t), W_R)
         DO i_cell_R = 1, bs_env%nimages_scf_desymm
   
            cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
   
            DO i_cell_R2 = 1, bs_env%nimages_scf_desymm
               cell_R2 = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R2, 1:3)
               IF (cell_R2(1) == -cell_R(1) .AND. cell_R2(2) == -cell_R(2)) THEN
                  i_cell_R_inv = i_cell_R2
                  cell_R_inv = cell_R2
               END IF
            END DO
 
            DO P_RI = 1, MIN(bs_env%n_RI, 10)
               DO Q_RI = 1, MIN(bs_env%n_RI, 10)
                  IF (bs_env%para_env%mepos == 0) THEN
                     PRINT *, 'i_t =', i_t,'P =', P_RI, 'Q =', Q_RI, &
                       'R =', cell_R(1:2), 'W^tr_PQ^R before *M =', W_R(P_RI, Q_RI, i_cell_R)
                  END IF
               END DO
            END DO
 
            DO P_RI = 1, MIN(bs_env%n_RI, 10)
               DO Q_RI = 1, MIN(bs_env%n_RI, 10)
                  IF (bs_env%para_env%mepos == 0) THEN
                     PRINT *, 'i_t =', i_t,'P =', P_RI, 'Q =', Q_RI, &
                       'R =', cell_R(1:2), 'W_(P,Q) after mult with M^-1 =', &
                    MWM_R(P_RI, Q_RI, i_cell_R)
                     IF (ABS(MWM_R(P_RI, Q_RI, i_cell_R)) > 1E-6_dp .AND. &
                         ABS(MWM_R(P_RI, Q_RI, i_cell_R) - &
                             MWM_R(Q_RI, P_RI, i_cell_R_inv))/ &
                         ABS(MWM_R(P_RI, Q_RI, i_cell_R)) > 1E-2_dp) THEN
                        PRINT *, 'ERROR'
                        PRINT *, 'i_t =', i_t, 'P =', Q_RI, 'Q =', P_RI, &
                                 'R =', cell_R_inv(1:2), 'Ŵ_PQ^R =', &
                           MWM_R(Q_RI, P_RI, i_cell_R_inv)
                     END IF
                  END IF
               END DO
            END DO
   
        END DO ! i_cell_R

      END DO ! i_t

! END JW 2del

      CALL timestop(handle)

   END SUBROUTINE compute_W_real_space

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param M_inv_V_sqrt ...
!> \param M_inv ...
!> \param V_sqrt ...
! **************************************************************************************************
   SUBROUTINE compute_Minv_and_Vsqrt(bs_env, qs_env, M_inv_V_sqrt, M_inv, V_sqrt, V_tr_BvK_sqrt, &
                                     M_inv_V_tr_BvK_sqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: M_inv_V_sqrt, M_inv, V_sqrt, &
      M_inv_V_tr_BvK_sqrt, V_tr_BvK_sqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Minv_and_Vsqrt'

      INTEGER                                            :: handle, ikp, ikp_local, n_RI, nkp, &
                                                            nkp_local, re_im, i_P, j_Q, i_cell, &
                                                            nkp_orig
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: M_tmp, V_tmp, M_R
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_V_kp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      nkp = bs_env%nkp_chi_eps_W_orig_plus_extra
      nkp_orig = bs_env%nkp_chi_eps_W_orig
      n_RI = bs_env%n_RI

      nkp_local = 0
      DO ikp = 1, nkp
         ! trivial parallelization over k-points
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE
         nkp_local = nkp_local + 1
      END DO

      ALLOCATE (M_inv_V_sqrt(n_RI, n_RI, nkp_local), M_inv(n_RI, n_RI, nkp_local), &
                V_sqrt(n_RI, n_RI, nkp_local), V_tr_BvK_sqrt(n_RI, n_RI, nkp_local), &
                M_inv_V_tr_BvK_sqrt(n_RI, n_RI, nkp_local))
      ALLOCATE (V_tmp(n_RI, n_RI, 1), M_tmp(n_RI, n_RI, 1))
      M_inv_V_sqrt(:, :, :) = z_zero
      M_inv(:, :, :) = z_zero
      V_sqrt(:, :, :) = z_zero
      V_tr_BvK_sqrt(:, :, :) = z_zero

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      NULLIFY (mat_V_kp)
      ALLOCATE (mat_V_kp(1:nkp, 2))
      DO ikp = 1, nkp
         DO re_im = 1, 2
            NULLIFY (mat_V_kp(ikp, re_im)%matrix)
            ALLOCATE (mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_create(mat_V_kp(ikp, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
            CALL dbcsr_reserve_all_blocks(mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_set(mat_V_kp(ikp, re_im)%matrix, 0.0_dp)
         END DO ! re_im
      END DO ! ikp

     ! 1. 2c Coulomb integrals for the first "original" k-point grid
      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
      CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, basis_type="RI_AUX", &
                                      cell=cell, particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=1, ikp_end=nkp_orig)

     ! 2. 2c Coulomb integrals for the second "extrapolation" k-point grid
      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_extra
      CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, basis_type="RI_AUX", &
                                      cell=cell, particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=nkp_orig+1, ikp_end=nkp)

      ALLOCATE (fm_V_kp(1))
      CALL cp_fm_create(fm_V_kp(1), bs_env%fm_RI_RI%matrix_struct)

      ikp_local = 0

      DO ikp = 1, nkp
         DO re_im = 1, 2

            CALL copy_dbcsr_to_fm(mat_V_kp(ikp, re_im)%matrix, fm_V_kp(1))
            CALL dbcsr_deallocate_matrix(mat_V_kp(ikp, re_im)%matrix)

            CALL bs_env%para_env%sync()

            CALL fm_to_local_array(fm_V_kp, V_tmp)

            IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

            SELECT CASE (re_im)
            CASE (1)
               ikp_local = ikp_local + 1
               V_sqrt(:, :, ikp_local) = V_sqrt(:, :, ikp_local) + z_one*V_tmp(:, :, 1)
            CASE (2)
               V_sqrt(:, :, ikp_local) = V_sqrt(:, :, ikp_local) + gaussi*V_tmp(:, :, 1)
               ! compute sqrt of V_PQ(k)
               CALL power(V_sqrt(:, :, ikp_local), 0.5_dp, 0.0_dp)
            END SELECT

         END DO
      END DO

      CALL cp_fm_release(fm_V_kp(1))
      DEALLOCATE (fm_V_kp)

      ! now get M^-1(k) and M^-1(k)*V^0.5(k)

      ! compute M^R_PQ = <phi_P,0|V^tr(rc=3Å)|phi_Q,R> for RI metric
      CALL get_V_tr_R(M_R, bs_env%ri_metric, bs_env%regularization_RI, bs_env, qs_env)

      ! JW 2del: 
      DO i_cell = 1, SIZE(M_R,3)
         PRINT *, '1. M_11^R for cell', i_cell, 'M =', M_R(1,1,i_cell)
      END DO

      ikp_local = 0
      DO ikp = 1, nkp

         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

         ikp_local = ikp_local + 1

         ! M(k) = sum_R e^ikR M^R
         CALL trafo_rs_to_ikp(M_R, M_inv(:, :, ikp_local), &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))

         ! JW 2del
         DO i_P = 1, MIN(5,bs_env%n_RI)
            DO j_Q = 1, MIN(5,bs_env%n_RI)
              PRINT *, 'for Minv*Vsqrt: M_PQ at xkp = ', &
                       bs_env%kpoints_chi_eps_W%xkp(1:3,ikp_local), &
                       i_P, j_Q, M_inv(i_P, j_Q, ikp_local)
              END DO
         END DO

         ! invert M_PQ(k)
         CALL power(M_inv(:, :, ikp_local), -1.0_dp, 0.0_dp)

         ! JW 2del
         DO i_P = 1, MIN(5,bs_env%n_RI)
            DO j_Q = 1, MIN(5,bs_env%n_RI)
              PRINT *, 'for Minv*Vsqrt: M_inv_PQ at xkp = ', &
                       bs_env%kpoints_chi_eps_W%xkp(1:3,ikp_local), &
                       i_P, j_Q, M_inv(i_P, j_Q, ikp_local)
              END DO
         END DO

         ! compute M^-1(k)*V^0.5(k)
         CALL ZGEMM("N", "C", n_RI, n_RI, n_RI, z_one, M_inv(:, :, ikp_local), n_RI, &
                    V_sqrt(:, :, ikp_local), n_RI, z_zero, M_inv_V_sqrt(:, :, ikp_local), n_RI)

      END DO ! ikp

      ! now get V^tr,BvK^0.5(k)

      ! compute V^tr,BvK^R_PQ = <phi_P,0|V^tr(rc=BvK cell size / 2)|phi_Q,R> stored in M_R
      CALL get_V_tr_R(M_R, bs_env%trunc_coulomb, 0.0_dp, bs_env, qs_env)

      ikp_local = 0
      DO ikp = 1, nkp

         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

         ikp_local = ikp_local + 1

         ! V^tr,BvK(k) = sum_R e^ikR V^tr,BvK^R
         CALL trafo_rs_to_ikp(M_R, V_tr_BvK_sqrt(:, :, ikp_local), &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))

         ! square root to obtain V^tr,BvK^0.5(k)
         CALL power(V_tr_BvK_sqrt(:, :, ikp_local), 0.5_dp, 0.0_dp)

         ! compute M^-1(k)*V^0.5(k)
         CALL ZGEMM("N", "C", n_RI, n_RI, n_RI, z_one, M_inv(:, :, ikp_local), n_RI, &
                    V_tr_BvK_sqrt(:, :, ikp_local), n_RI, z_zero, &
                    M_inv_V_tr_BvK_sqrt(:, :, ikp_local), n_RI)

      END DO ! ikp

      CALL timestop(handle)

   END SUBROUTINE compute_Minv_and_Vsqrt

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param exponent ...
!> \param eps ...
! **************************************************************************************************
   SUBROUTINE power(matrix, exponent, eps)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: matrix
      REAL(KIND=dp)                                      :: exponent, eps

      CHARACTER(len=*), PARAMETER                        :: routineN = 'power'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenvectors
      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: work
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: A
      INTEGER                                            :: handle, i, info, liwork, lrwork, lwork, n
      INTEGER, DIMENSION(:), POINTER                     :: iwork
      REAL(KIND=dp)                                      :: pos_eval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rwork

      CALL timeset(routineN, handle)

      ! code by Ole Schütt
      IF (SIZE(matrix, 1) /= SIZE(matrix, 2)) CPABORT("expected square matrix")

      ! make matrix perfectly Hermitian
      matrix(:, :) = 0.5_dp*(matrix(:, :) + CONJG(TRANSPOSE(matrix(:, :))))

      n = SIZE(matrix, 1)
      ALLOCATE (iwork(1), rwork(1), work(1), A(n, n), eigenvalues(n), eigenvectors(n, n))

      A(:, :) = matrix ! ZHEEVD will overwrite A
      ! work space query
      lwork = -1
      lrwork = -1
      liwork = -1

      CALL ZHEEVD('V', 'U', n, A(1, 1), n, eigenvalues(1), &
                  work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
      lwork = INT(REAL(work(1), dp))
      lrwork = INT(REAL(rwork(1), dp))
      liwork = iwork(1)

      DEALLOCATE (iwork, rwork, work)
      ALLOCATE (iwork(liwork))
      iwork(:) = 0
      ALLOCATE (rwork(lrwork))
      rwork(:) = 0.0_dp
      ALLOCATE (work(lwork))
      work(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      CALL ZHEEVD('V', 'U', n, A(1, 1), n, eigenvalues(1), &
                  work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)

      eigenvectors = A

      IF (info /= 0) CPABORT("diagonalization failed")

      DO i = 1, n
         IF (eigenvalues(i) > eps) THEN
            pos_eval = (eigenvalues(i))**(0.5_dp*exponent)
         ELSE
            pos_eval = z_zero
         END IF
         eigenvectors(:, i) = eigenvectors(:, i)*pos_eval
      END DO

      CALL ZGEMM("N", "C", n, n, n, z_one, eigenvectors, n, eigenvectors, n, z_zero, matrix, n)

      DEALLOCATE (iwork, rwork, work, A, eigenvalues, eigenvectors)

      CALL timestop(handle)

   END SUBROUTINE power

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE add_on_diag(matrix, alpha)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: matrix
      COMPLEX(KIND=dp)                                   :: alpha

      CHARACTER(len=*), PARAMETER                        :: routineN = 'add_on_diag'

      INTEGER                                            :: handle, i, n

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(n == SIZE(matrix, 2))

      DO i = 1, n
         matrix(i, i) = matrix(i, i) + alpha
      END DO

      CALL timestop(handle)

   END SUBROUTINE add_on_diag

   SUBROUTINE truncate_W_and_mult_W_with_Minv(W_R, MWM_R, bs_env, qs_env, j_w)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'truncate_W_and_mult_W_with_Minv'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chi_k_w, eps_k_w, W_k, M_inv, work, M
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: M_inv_V_sqrt, V_sqrt, &
                 V_tr_BvK_sqrt, M_inv_V_tr_BvK_sqrt
      INTEGER                                            :: handle, i_cell, i_cell_R, i_cell_R2, &
                                                            i_cell_R_inv, i_t, ikp, ikp_local, &
                                                            j_w, n_RI, nkp, P_RI, Q_RI, i_P, j_Q, &
                                                            nkp_local, P, Q
      INTEGER, DIMENSION(3)                              :: cell_R, cell_R2, cell_R_inv
      REAL(KIND=dp)                                      :: cell(3), el, freq_j, kR, time_i, &
                                                            weight_ij, wkp, xkp(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: chi_R, W_R, M_R, MWM_R

      CALL timeset(routineN, handle)

      ! 1.  set W_PQ^R(iω) to zero if distance of φ_P (in cell 0) and φ_Q (in cell R) is 
      !     larger than half the Born-von-Karman supercell (as usually done in periodic HFX);
      !     we call result W_PQ^R^trunc(iω)
! JW 2del NO TRUNCATION AND REMOVE TRUNCATE FROM SUBROUTINE NAMES!!!  
! JW 2del    CALL truncate_W(W_R, bs_env, qs_env, j_w)

      ! 2.  compute Ŵ(iω,k) = M^-1(k) W^trunc(iω,k) M^-1(k) from W^R^trunc(iω) and trafo to Ŵ^R(iω)

      ! 2a. compute M^R again
      CALL get_V_tr_R(M_R, bs_env%ri_metric, bs_env%regularization_RI, bs_env, qs_env)

      ! JW 2del:
      DO i_cell = 1, SIZE(M_R,3)
         PRINT *, '2. M_11^R for cell', i_cell, 'M =', M_R(1,1,i_cell)
      END DO

      n_RI = bs_env%n_RI
! JW 2del: M
      ALLOCATE(M_inv(n_RI, n_RI), W_k(n_RI, n_RI), work(n_RI, n_RI), M(n_RI, n_RI))
      MWM_R(:,:,:) = 0.0_dp
      DO ikp = 1, bs_env%nkp_scf_desymm

         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

         ! 2b. M(k) = sum_R e^ikR M^R
         CALL trafo_rs_to_ikp(M_R, M_inv, &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))

         PRINT *, 'xkp =', bs_env%kpoints_scf_desymm%xkp(1:2, ikp), &
                  'M_11 =', M_inv(1,1)

         IF(ikp == 1) THEN

            DO P = 1, bs_env%n_RI 
            DO Q = 1, bs_env%n_RI
               PRINT *, 'M_PQ: P =', P, Q, 'M_PQ =', REAL(M_inv(P,Q))
            END DO
            END DO

         END IF

         ! JW 2del
         M(:,:) = M_inv(:,:)

         ! 2c. invert M_PQ(k)
         CALL power(M_inv, -1.0_dp, 0.0_dp)

         PRINT *, 'xkp =', bs_env%kpoints_scf_desymm%xkp(1:2, ikp), &
                  'M_inv_11 =', M_inv(1,1), 'W_11 =', W_k(1,1)

         ! JW 2del
         IF(ikp == 1) THEN

            CALL ZGEMM("N", "N", n_RI, n_RI, n_RI, z_one, M_inv, n_RI, M, n_RI, z_zero, work, n_RI)

            DO P = 1, bs_env%n_RI
            DO Q = 1, bs_env%n_RI
               PRINT *, 'M_inv_PQ: P =', P, Q, 'M_PQ =', M_inv(P,Q), 'M*M_inv_PQ =', work(P,Q)
            END DO
            END DO

         END IF

         ! 2d. W^trunc(k) = sum_R e^ikR W^trunc^R
         CALL trafo_rs_to_ikp(W_R, W_k, &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))

         ! 2e. M^-1(k) W^trunc(k)
         CALL ZGEMM("N", "N", n_RI, n_RI, n_RI, z_one, M_inv, n_RI, W_k, n_RI, z_zero, work, n_RI)

         ! 2f. Ŵ(k) = M^-1(k) W^trunc(k) M^-1(k)
         CALL ZGEMM("N", "N", n_RI, n_RI, n_RI, z_one, work, n_RI, M_inv, n_RI, z_zero, W_k, n_RI)

         ! 2g. Ŵ^R = sum_k w_k e^(-ikR) Ŵ^(k)
         CALL add_ikp_to_all_rs(W_k, MWM_R, bs_env%kpoints_scf_desymm, ikp)

      END DO ! ikp

      CALL bs_env%para_env%sync()
      CALL bs_env%para_env%sum(MWM_R)

      CALL timestop(handle)

   END SUBROUTINE truncate_W_and_mult_W_with_Minv

   SUBROUTINE  truncate_W(W_R, bs_env, qs_env, j_w)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'truncate_W'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chi_k_w, eps_k_w, W_k, M_inv, work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: M_inv_V_sqrt, V_sqrt, &
                 V_tr_BvK_sqrt, M_inv_V_tr_BvK_sqrt
      INTEGER                                            :: handle, i_cell, i_cell_R, i_cell_R2, &
                                                            i_cell_R_inv, i_t, ikp, ikp_local, &
                                                            j_w, n_RI, nkp, P_RI, Q_RI, i_P, j_Q, &
                                                            nkp_local, i_counter, P, Q, &
     P_atom, Q_atom, n_cells
      INTEGER, DIMENSION(3)                              :: cell_R, cell_R2, cell_R_inv, cell_coord
      REAL(KIND=dp)                                      :: el, freq_j, kR, time_i, d_PQ, r_PQ(3), &
                                                            weight_ij, wkp, xkp(3), cell_vector(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: chi_R, W_R, M_R, MWM_R, W_R_tmp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf_index

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      n_cells = bs_env%nimages_scf_desymm

      ALLOCATE (atom_from_bf_index(n_RI))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf_index, n_RI, "RI_AUX")

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      ALLOCATE(W_R_tmp(n_RI, n_RI, n_cells))
      W_R_tmp(:,:,:) = 0.0_dp

      i_counter = 0
      DO P = 1, n_RI
        DO Q = 1, n_RI 
          DO i_cell = 1, n_cells

            ! trivial parallelization
            i_counter = i_counter + 1
            IF (MODULO(i_counter, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE 

            P_atom = atom_from_bf_index(P)
            Q_atom = atom_from_bf_index(Q)

            cell_coord(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell, 1:3)

            cell_vector(1:3) = MATMUL(hmat, REAL(cell_coord, dp))

            r_PQ(1:3) = pbc(particle_set(P_atom)%r(1:3), cell) - &
                        pbc(particle_set(Q_atom)%r(1:3), cell) + cell_vector(1:3)

            d_PQ = SQRT(r_PQ(1)**2 + r_PQ(2)**2 + r_PQ(3)**2)

            ! truncate the screened Coulomb interaction after half the BvK cell radius
            IF(d_PQ < 0.5_dp*bs_env%trunc_coulomb%cutoff_radius) THEN
              W_R_tmp(P, Q, i_cell) = W_R(P, Q, i_cell)
            END IF

          END DO
        END DO
      END DO

      CALL bs_env%para_env%sync()
      CALL bs_env%para_env%sum(W_R_tmp)

      DO P = 1, n_RI
        DO Q = 1, n_RI
          DO i_cell = 1, n_cells

            cell_coord = bs_env%kpoints_scf_desymm%index_to_cell(i_cell, 1:3)

            IF(bs_env%para_env%mepos == 0 .AND. P == 1 .AND. Q == 1 .AND. j_w == 1) THEN
               PRINT *, 'P, Q =', P, Q, 'R =', cell_coord, 'Wtr =', W_R_tmp(P, Q, i_cell), &
                        'W =', W_R(P, Q, i_cell)
            END IF

          END DO
        END DO
      END DO

      W_R(:,:,:) = W_R_tmp(:,:,:)

      CALL timestop(handle)

   END SUBROUTINE truncate_W

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, i_ao_1, i_ao_2, i_cell_R, &
                                                            i_task_local, ispin, nimages_3c
      INTEGER(KIND=int_8)                                :: n_double_in_3c
      INTEGER, DIMENSION(3)                              :: cell_R
      LOGICAL                                            :: memory_fully_used
      REAL(KIND=dp)                                      :: t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: D_S, Mi_Vtr_Mi_R, Sigma_x_R
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_D, M_V
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c_D, store_int_3c_V

      CALL timeset(routineN, handle)

      t1 = m_walltime()
      n_double_in_3c = INT(0, KIND=int_8)
      nimages_3c = bs_env%nimages_3c
      ALLOCATE (M_D(nimages_3c, nimages_3c))
      ALLOCATE (M_V(nimages_3c, nimages_3c))
      ALLOCATE (store_int_3c_V(nimages_3c, nimages_3c, bs_env%n_atom))
      ALLOCATE (store_int_3c_D(nimages_3c, nimages_3c, bs_env%n_atom))

      ! density matrix D_µλ^S, S: cell; Σ^x_λσ^R, R: cell, replicated on each MPI process
      ALLOCATE (D_S(bs_env%n_ao, bs_env%n_ao, bs_env%nimages_scf_desymm))
      ALLOCATE (Sigma_x_R(bs_env%n_ao, bs_env%n_ao, bs_env%nimages_scf_desymm))
      Sigma_x_R(:, :, :) = 0.0_dp

      ! V^tr_PQ^R = <phi_P,0|V^tr|phi_Q,R>, V^tr(k) = sum_R e^ikR V^tr^R
      ! M(k) = sum_R e^ikR M^R, M(k) -> M^-1(k) -> Ṽ^tr(k) = M^-1(k) * V^tr(k) * M^-1(k)
      !                                         -> Ṽ^tr_PQ^R = sum_k w_k e^-ikR Ṽ^tr_PQ(k)
      CALL get_MinvVtrMinv_R(Mi_Vtr_Mi_R, bs_env, qs_env)

      memory_fully_used = .FALSE.

      ! Σ^x_λσ^R = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) D_µν^S2       ]
      !                       [ sum_QR2 (σR νS1    | QR1-R2) Ṽ^tr_PQ^R2 ]
      DO ispin = 1, bs_env%n_spin

         ! compute D^S(iτ) for cell S from D_µν(k) = sum_n^occ C^*_µn(k) C_νn(k):
         ! trafo k-point k -> cell S: D_µν^S = sum_k w_k D_µν(k) e^(ikS)
         CALL G_occ_vir(bs_env, 0.0_dp, D_S, ispin, occ=.TRUE., vir=.FALSE.)

         ! JW 2del
         IF (bs_env%para_env%mepos == 0) THEN
            DO i_cell_R = 1, bs_env%nimages_scf_desymm
               cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
               DO i_ao_1 = 1, bs_env%n_ao
                  DO i_ao_2 = 1, bs_env%n_ao
                     PRINT *, 'µ =', i_ao_1, 'ν =', i_ao_2, 'R =', cell_R(1:2), &
                        'D_µν^R =', D_S(i_ao_1, i_ao_2, i_cell_R)
                  END DO
               END DO
            END DO
         END IF

         IF (bs_env%para_env%mepos == 0) THEN
            DO i_cell_R = 1, bs_env%nimages_scf_desymm
               cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
               DO i_ao_1 = 1, MIN(bs_env%n_RI, 20)
                  DO i_ao_2 = 1, MIN(bs_env%n_RI, 20)
                     PRINT *, 'P =', i_ao_1, 'Q =', i_ao_2, 'R =', cell_R(1:2), &
                        'Mi_Vtr_Mi_PQ^R =', Mi_Vtr_Mi_R(i_ao_1, i_ao_2, i_cell_R)
                  END DO
               END DO
            END DO
         END IF
         ! end 2del

         DO i_task_local = 1, bs_env%n_tasks_local

            ! M^V_σ0,νS1,PR1 = sum_QR2 ( σ0 νS1 | QR1-R2 ) Ṽ^tr_QP^R2 for i_task_local
            CALL compute_3c_and_contract_W(M_V, Mi_Vtr_Mi_R, store_int_3c_V, i_task_local, &
                                           n_double_in_3c, memory_fully_used, bs_env, qs_env)
            ! M^D_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) D_µν^S2
            ! Σ_λσ^R = sum_PR1νS1 M^D_λ0,νS1,PR1 * M^V_σR,νS1,PR1 for i_task_local, where
            !                                      M^V_σR,νS1,PR1 = M^V_σ0,νS1-R,PR1-R
            CALL contract_to_Sigma(Sigma_x_R, D_S, M_V, store_int_3c_D, i_task_local, &
                                   n_double_in_3c, memory_fully_used, bs_env, qs_env, &
                                   occ=.TRUE., vir=.FALSE.)

         END DO ! i_task_local: atom_P, atom_ν, i_cell_Delta_R

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(Sigma_x_R)

         CALL local_array_to_fm(Sigma_x_R, bs_env%fm_Sigma_x_R(:))

         ! JW 2del
         IF (bs_env%para_env%mepos == 0) THEN
            DO i_cell_R = 1, bs_env%nimages_scf_desymm
               cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
               DO i_ao_1 = 1, bs_env%n_ao
                  DO i_ao_2 = 1, bs_env%n_ao
                     PRINT *, 'λ =', i_ao_1, 'σ =', i_ao_2, 'R =', cell_R(1:2), &
                        'Σ_λσ^R =', Sigma_x_R(i_ao_1, i_ao_2, i_cell_R)
                  END DO
               END DO
            END DO
         END IF
         ! end 2del

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief ...
!> \param Mi_Vtr_Mi_R ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_MinvVtrMinv_R(Mi_Vtr_Mi_R, bs_env, qs_env)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Mi_Vtr_Mi_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_MinvVtrMinv_R'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: M_inv_V_tr_kp, M_kp, Mi_Vtr_Mi_kp, &
                                                            V_tr_kp
      INTEGER                                            :: handle, ikp, img, n_RI, nimages_scf, &
                                                            nkp_scf, P, Q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: M_R, V_tr_R

      CALL timeset(routineN, handle)

      nimages_scf = bs_env%nimages_scf_desymm
      nkp_scf = bs_env%kpoints_scf_desymm%nkp
      n_RI = bs_env%n_RI

      CALL get_V_tr_R(V_tr_R, bs_env%trunc_coulomb, 0.0_dp, bs_env, qs_env)
      CALL get_V_tr_R(M_R, bs_env%ri_metric, bs_env%regularization_RI, bs_env, qs_env)

      ALLOCATE (V_tr_kp(n_RI, n_RI), M_kp(n_RI, n_RI), M_inv_V_tr_kp(n_RI, n_RI), &
                Mi_Vtr_Mi_kp(n_RI, n_RI), Mi_Vtr_Mi_R(n_RI, n_RI, nimages_scf))
      Mi_Vtr_Mi_R(:, :, :) = 0.0_dp

      DO ikp = 1, nkp_scf
         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE
         ! V_tr(k) = sum_R e^ikR V_tr^R
         CALL trafo_rs_to_ikp(V_tr_R, V_tr_kp, bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))
         ! M(k)    = sum_R e^ikR M^R
         CALL trafo_rs_to_ikp(M_R, M_kp, bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))
         ! M(k) -> M^-1(k)
         CALL power(M_kp, -1.0_dp, 0.0_dp)
         ! M^-1(k) * V_tr(k)
         CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, M_kp, n_RI, &
                    V_tr_kp, n_RI, z_zero, M_inv_V_tr_kp, n_RI)
         ! Ṽ(k) = M^-1(k) * V_tr(k) * M^-1(k)
         CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, M_inv_V_tr_kp, n_RI, &
                    M_kp, n_RI, z_zero, Mi_Vtr_Mi_kp, n_RI)
         ! Ṽ^R = sum_k w_k e^-ikR Ṽ(k)
         CALL add_ikp_to_all_rs(Mi_Vtr_Mi_kp, Mi_Vtr_Mi_R, bs_env%kpoints_scf_desymm, ikp)
      END DO
      CALL bs_env%para_env%sync()
      CALL bs_env%para_env%sum(Mi_Vtr_Mi_R)

      ! JW 2del
      IF (bs_env%para_env%num_pe == 1) THEN
         DO img = 1, nimages_scf
            DO P = 1, bs_env%n_RI
               DO Q = 1, bs_env%n_RI
                  PRINT *, 'R =', bs_env%kpoints_scf_desymm%index_to_cell(img, 1:2), &
                     'P, Q =', P, Q, 'V_tr, MiV_trMi =', V_tr_R(P, Q, img), Mi_Vtr_Mi_R(P, Q, img)
               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_MinvVtrMinv_R

! **************************************************************************************************
!> \brief ...
!> \param V_tr_R ...
!> \param pot_type ...
!> \param regularization_RI ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_V_tr_R(V_tr_R, pot_type, regularization_RI, bs_env, qs_env)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: V_tr_R
      TYPE(libint_potential_type)                        :: pot_type
      REAL(KIND=dp)                                      :: regularization_RI
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_V_tr_R'

      INTEGER                                            :: handle, img, nimages_scf_desymm
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_V_tr_R
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_V_tr_R
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (sab_RI, dist_2d)

      CALL get_qs_env(qs_env=qs_env, &
                      blacs_env=blacs_env, &
                      distribution_2d=dist_2d, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      ALLOCATE (sizes_RI(bs_env%n_atom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=bs_env%basis_set_RI)
      CALL build_2c_neighbor_lists(sab_RI, bs_env%basis_set_RI, bs_env%basis_set_RI, &
                                   pot_type, "2c_nl_RI", qs_env, sym_ij=.FALSE., &
                                   dist_2d=dist_2d)
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(sizes_RI)))
      ALLOCATE (col_bsize(SIZE(sizes_RI)))
      row_bsize(:) = sizes_RI
      col_bsize(:) = sizes_RI

      nimages_scf_desymm = bs_env%nimages_scf_desymm
      ALLOCATE (mat_V_tr_R(nimages_scf_desymm))
      CALL dbcsr_create(mat_V_tr_R(1), "(RI|RI)", dbcsr_dist, dbcsr_type_no_symmetry, &
                        row_bsize, col_bsize)
      DEALLOCATE (row_bsize, col_bsize)

      DO img = 2, nimages_scf_desymm
         CALL dbcsr_create(mat_V_tr_R(img), template=mat_V_tr_R(1))
      END DO

      CALL build_2c_integrals(mat_V_tr_R, 0.0_dp, qs_env, sab_RI, bs_env%basis_set_RI, &
                              bs_env%basis_set_RI, pot_type, do_kpoints=.TRUE., &
                              ext_kpoints=bs_env%kpoints_scf_desymm, &
                              regularization_RI=regularization_RI)

      ALLOCATE (fm_V_tr_R(nimages_scf_desymm))
      DO img = 1, nimages_scf_desymm
         CALL cp_fm_create(fm_V_tr_R(img), bs_env%fm_RI_RI%matrix_struct)
         CALL copy_dbcsr_to_fm(mat_V_tr_R(img), fm_V_tr_R(img))
         CALL dbcsr_release(mat_V_tr_R(img))
      END DO

      IF(.NOT. ALLOCATED(V_tr_R)) THEN
         ALLOCATE (V_tr_R(bs_env%n_RI, bs_env%n_RI, nimages_scf_desymm))
      END IF

      CALL fm_to_local_array(fm_V_tr_R, V_tr_R)

      CALL cp_fm_release(fm_V_tr_R)
      CALL dbcsr_distribution_release(dbcsr_dist)
      CALL release_neighbor_list_sets(sab_RI)

      CALL timestop(handle)

   END SUBROUTINE get_V_tr_R

! **************************************************************************************************
!> \brief ...
!> \param M_W ...
!> \param W ...
!> \param store_int_3c ...
!> \param i_task_local ...
!> \param n_double_in_3c ...
!> \param memory_fully_used ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_3c_and_contract_W(M_W, W, store_int_3c, i_task_local, &
                                        n_double_in_3c, memory_fully_used, bs_env, qs_env)
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_W
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: W
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c
      INTEGER                                            :: i_task_local
      INTEGER(KIND=int_8)                                :: n_double_in_3c
      LOGICAL                                            :: memory_fully_used
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_and_contract_W'

      INTEGER :: atom_nu, atom_P, handle, handle2, i_ao_end, i_ao_size, i_ao_start, &
         i_cell_Delta_R, i_cell_j, i_cell_k, i_cell_R1, i_cell_R2, i_cell_S1, i_P, i_Q, i_RI_end, &
         i_RI_size, i_RI_start, img, n_ao, n_RI
      INTEGER, DIMENSION(3)                              :: cell_R1, cell_R1_minus_cell_R2, cell_R2, &
                                                            cell_S1
      LOGICAL                                            :: cell_found, is_S1_in_3c_cells
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c, M_tmp

      CALL timeset(routineN, handle)

      atom_nu = bs_env%task_atom_1_atom_2_img(i_task_local, 1)
      atom_P = bs_env%task_atom_1_atom_2_img(i_task_local, 2)
      i_cell_Delta_R = bs_env%task_atom_1_atom_2_img(i_task_local, 3)

      i_ao_start = bs_env%i_ao_start_from_atom(atom_nu)
      i_ao_end = bs_env%i_ao_end_from_atom(atom_nu)
      i_ao_size = i_ao_end - i_ao_start + 1

      i_RI_start = bs_env%i_RI_start_from_atom(atom_P)
      i_RI_end = bs_env%i_RI_end_from_atom(atom_P)
      i_RI_size = i_RI_end - i_RI_start + 1

      n_ao = bs_env%n_ao
      n_RI = bs_env%n_RI

      ALLOCATE (int_3c(1:n_ao, 1:i_ao_size, 1:n_RI))
      ALLOCATE (M_tmp(1:n_ao, 1:i_ao_size, 1:i_RI_size))

      DO i_cell_R1 = 1, bs_env%nimages_3c

         CALL get_second_R(i_cell_Delta_R, i_cell_R1, bs_env, is_S1_in_3c_cells, i_cell_S1)

         IF (.NOT. is_S1_in_3c_cells) CYCLE

         ! deallocate if from previous i_task_local with another atom_P, atom_nu
         IF (ALLOCATED(M_W(i_cell_S1, i_cell_R1)%data_3)) THEN
            DEALLOCATE (M_W(i_cell_S1, i_cell_R1)%data_3)
         END IF
         ALLOCATE (M_W(i_cell_S1, i_cell_R1)%data_3(1:bs_env%n_ao, 1:i_ao_size, 1:i_RI_size))
         M_W(i_cell_S1, i_cell_R1)%data_3(:, :, :) = 0.0_dp

         M_tmp(:, :, :) = 0.0_dp

         DO i_cell_R2 = 1, bs_env%nimages_scf_desymm

            cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
            cell_S1(1:3) = bs_env%index_to_cell_3c(i_cell_S1, 1:3)
            cell_R2(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R2, 1:3)
            cell_R1_minus_cell_R2(1:3) = cell_R1(1:3) - cell_R2(1:3)

            CALL is_cell_in_index_to_cell(cell_R1_minus_cell_R2, bs_env%index_to_cell_3c, cell_found)

            IF (.NOT. cell_found) CYCLE

            i_cell_k = bs_env%cell_to_index_3c(cell_S1(1), cell_S1(2), cell_S1(3))
            i_cell_j = bs_env%cell_to_index_3c(cell_R1_minus_cell_R2(1), cell_R1_minus_cell_R2(2), &
                                               cell_R1_minus_cell_R2(3))

            IF(bs_env%nblocks_3c(i_cell_j, i_cell_k) == 0) CYCLE

            IF (ALLOCATED(store_int_3c(i_cell_j, i_cell_k, atom_nu)%data_3)) THEN
               int_3c(:, :, :) = store_int_3c(i_cell_j, i_cell_k, atom_nu)%data_3(:, :, :)
            ELSE

               ! compute 3-c int. ( σ (all at., cell 0) , ν (atom_ν, cell S1) | Q (all at., R1-R2) )
               ! ("|": truncated Coulomb operator), inside build_3c_integral_block: (j k | i)
               CALL build_3c_integral_block(int_3c, qs_env, bs_env%ri_metric, &
                                       basis_j=bs_env%basis_set_AO, &
                                       basis_k=bs_env%basis_set_AO, &
                                       basis_i=bs_env%basis_set_RI, &
                                       cell_k=cell_S1, &
                                       cell_i=cell_R1_minus_cell_R2, &
                                       atom_k=atom_nu, &
                                       j_bf_start_from_atom=bs_env%i_ao_start_from_atom, &
                                       i_bf_start_from_atom=bs_env%i_RI_start_from_atom)

               IF (n_double_in_3c < bs_env%avail_doubles_per_proc) THEN
                  ALLOCATE (store_int_3c(i_cell_j, i_cell_k, atom_nu)% &
                            data_3(1:n_ao, 1:i_ao_size, 1:n_RI))
                  store_int_3c(i_cell_j, i_cell_k, atom_nu)%data_3(:, :, :) = int_3c(:, :, :)
                  n_double_in_3c = n_double_in_3c + &
                                   INT(i_ao_size*n_ao, KIND=int_8)*INT(n_RI, KIND=int_8)
               ELSE
                  ! print warning when memory is full for first time
                  IF(.NOT. memory_fully_used) CALL warning(bs_env)
                  memory_fully_used = .TRUE.
               END IF
            END IF

            CALL timeset(routineN//"_3c_x_W", handle2)

            ! M^W_σ0,νS1,PR1 = sum_QR2 ( σ0 νS1 | QR1-R2 ) W_QP^R2 for i_task_local
            DO i_P = 1, i_RI_size
               DO i_Q = 1, n_RI
                  M_tmp(:, :, i_P) = M_tmp(:, :, i_P) + &
                                     int_3c(:, :, i_Q)* &
                                     W(i_Q, i_P + i_RI_start - 1, i_cell_R2)
               END DO
            END DO

            CALL timestop(handle2)

         END DO ! i_cell_S2

         M_W(i_cell_S1, i_cell_R1)%data_3(:, :, :) = M_W(i_cell_S1, i_cell_R1)%data_3(:, :, :) + &
                                                     M_tmp(:, :, :)

      END DO ! i_cell_R1

      CALL timestop(handle)

   END SUBROUTINE compute_3c_and_contract_W

! **************************************************************************************************
!> \brief ...
!> \param Sigma_R ...
!> \param G_S ...
!> \param M_W ...
!> \param store_int_3c ...
!> \param i_task_local ...
!> \param n_double_in_3c ...
!> \param memory_fully_used ...
!> \param bs_env ...
!> \param qs_env ...
!> \param occ ...
!> \param vir ...
! **************************************************************************************************
   SUBROUTINE contract_to_Sigma(Sigma_R, G_S, M_W, store_int_3c, i_task_local, n_double_in_3c, &
                                memory_fully_used, bs_env, qs_env, occ, vir)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_R, G_S
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_W
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c
      INTEGER                                            :: i_task_local
      INTEGER(KIND=int_8)                                :: n_double_in_3c
      LOGICAL                                            :: memory_fully_used
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: occ, vir

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_to_Sigma'

      INTEGER :: atom_nu, atom_P, handle, handle2, i_ao_1, i_ao_2, i_ao_end, i_ao_size, &
         i_ao_start, i_cell_Delta_R, i_cell_i, i_cell_k, i_cell_R, i_cell_R1, i_cell_R1_minus_R, &
         i_cell_S1, i_cell_S1_minus_R, i_cell_S2, i_mu, i_nu, i_RI_end, i_RI_size, i_RI_start, &
         img, n_ao, n_RI, i_cell_j_check, i_cell_k_check
      INTEGER, DIMENSION(3)                              :: cell_R, cell_R1, cell_S1, &
                                                            cell_S1_plus_cell_S2, cell_S2, &
        cell_S1_p_S2_m_R1
      LOGICAL                                            :: cell_found, is_R1_minus_R_in_3c_cells, &
                                                            is_S1_in_3c_cells, &
                                                            is_S1_minus_R_in_3c_cells
      REAL(KIND=dp)                                      :: sign_Sigma, sum_GW
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: int_3c, M_G

      CALL timeset(routineN, handle)

      CPASSERT(occ .EQV. (.NOT. vir))
      IF (occ) sign_Sigma = -1.0_dp
      IF (vir) sign_Sigma = 1.0_dp

      atom_nu = bs_env%task_atom_1_atom_2_img(i_task_local, 1)
      atom_P = bs_env%task_atom_1_atom_2_img(i_task_local, 2)
      i_cell_Delta_R = bs_env%task_atom_1_atom_2_img(i_task_local, 3)

      i_ao_start = bs_env%i_ao_start_from_atom(atom_nu)
      i_ao_end = bs_env%i_ao_end_from_atom(atom_nu)
      i_ao_size = i_ao_end - i_ao_start + 1

      i_RI_start = bs_env%i_RI_start_from_atom(atom_P)
      i_RI_end = bs_env%i_RI_end_from_atom(atom_P)
      i_RI_size = i_RI_end - i_RI_start + 1

      n_ao = bs_env%n_ao
      n_RI = bs_env%n_RI

      ALLOCATE (int_3c(1:n_ao, 1:n_ao, 1:i_RI_size))
      ALLOCATE (M_G(1:n_ao, 1:i_ao_size, 1:i_RI_size))

      DO i_cell_R1 = 1, bs_env%nimages_3c

         CALL get_second_R(i_cell_Delta_R, i_cell_R1, bs_env, is_S1_in_3c_cells, i_cell_S1)

         IF (.NOT. is_S1_in_3c_cells) CYCLE

         M_G(:, :, :) = 0.0_dp

         DO i_cell_S2 = 1, bs_env%nimages_scf_desymm

            cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
            cell_S1(1:3) = bs_env%index_to_cell_3c(i_cell_S1, 1:3)
            cell_S2(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_S2, 1:3)
            cell_S1_plus_cell_S2(1:3) = cell_S1(1:3) + cell_S2(1:3)

! JW 2del
!            CALL is_cell_in_index_to_cell(cell_S1_plus_cell_S2, bs_env%index_to_cell_3c, cell_found)
!
!            IF (.NOT. cell_found) CYCLE

            cell_S1_p_S2_m_R1(1:3) = cell_S1(1:3) + cell_S2(1:3) - cell_R1(1:3)

            CALL is_cell_in_index_to_cell(cell_S1_p_S2_m_R1, bs_env%index_to_cell_3c, cell_found)

            IF (.NOT. cell_found) CYCLE

            i_cell_j_check = bs_env%cell_to_index_3c(-cell_R1(1), -cell_R1(2), -cell_R1(3))
            i_cell_k_check = bs_env%cell_to_index_3c(cell_S1_p_S2_m_R1(1), cell_S1_p_S2_m_R1(2), &
                                                     cell_S1_p_S2_m_R1(3))
            IF(bs_env%nblocks_3c(i_cell_j_check, i_cell_k_check) == 0) CYCLE

            i_cell_i = bs_env%cell_to_index_3c(cell_R1(1), cell_R1(2), cell_R1(3))
            i_cell_k = bs_env%cell_to_index_3c(cell_S1_plus_cell_S2(1), cell_S1_plus_cell_S2(2), &
                                               cell_S1_plus_cell_S2(3))

            IF (ALLOCATED(store_int_3c(i_cell_i, i_cell_k, atom_P)%data_3)) THEN
               int_3c(:, :, :) = store_int_3c(i_cell_i, i_cell_k, atom_P)%data_3(:, :, :)
            ELSE
               ! compute 3-c int. ( λ (all at., cell 0),  µ (all at., cell S1+S2) | P (atom_P, R1) )
               ! ("|": truncated Coulomb operator), inside build_3c_integral_block: (j k | i)
               CALL build_3c_integral_block(int_3c, qs_env, bs_env%ri_metric, &
                                       basis_j=bs_env%basis_set_AO, &
                                       basis_k=bs_env%basis_set_AO, &
                                       basis_i=bs_env%basis_set_RI, &
                                       cell_k=cell_S1_plus_cell_S2, &
                                       cell_i=cell_R1, &
                                       atom_i=atom_P, &
                                       j_bf_start_from_atom=bs_env%i_ao_start_from_atom, &
                                       k_bf_start_from_atom=bs_env%i_ao_start_from_atom)

               IF (n_double_in_3c < bs_env%avail_doubles_per_proc) THEN
                  ALLOCATE (store_int_3c(i_cell_i, i_cell_k, atom_P)% &
                            data_3(1:n_ao, 1:n_ao, 1:i_RI_size))
                  store_int_3c(i_cell_i, i_cell_k, atom_P)%data_3(:, :, :) = int_3c(:, :, :)
                  n_double_in_3c = n_double_in_3c + &
                                   INT(i_RI_size*n_ao, KIND=int_8)*INT(n_ao, KIND=int_8)
               ELSE
                  ! print warning when memory is full for first time
                  IF(.NOT. memory_fully_used) CALL warning(bs_env)
                  memory_fully_used = .TRUE.
               END IF
            END IF

            CALL timeset(routineN//"_3c_x_G", handle2)

            ! M_λ0,νS1,PR1 = sum_µS2 ( λ0 µS1-S2 | PR1 ) G^occ/vir_µν^S2(i|τ|) for i_task_local
            DO i_nu = 1, i_ao_size
               DO i_mu = 1, n_ao
                  M_G(:, i_nu, :) = M_G(:, i_nu, :) + &
                                    int_3c(:, i_mu, :)*G_S(i_nu + i_ao_start - 1, i_mu, i_cell_S2)
               END DO
            END DO

            CALL timestop(handle2)

         END DO ! i_cell_S2

         CALL timeset(routineN//"_contract", handle2)

         DO i_cell_R = 1, bs_env%nimages_scf_desymm

            cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)

            CALL get_Delta_R(i_cell_R1_minus_R, is_R1_minus_R_in_3c_cells, &
                             i_cell_R1, i_cell_R, bs_env)
            CALL get_Delta_R(i_cell_S1_minus_R, is_S1_minus_R_in_3c_cells, &
                             i_cell_S1, i_cell_R, bs_env)

            IF (.NOT. is_R1_minus_R_in_3c_cells) CYCLE
            IF (.NOT. is_S1_minus_R_in_3c_cells) CYCLE

            DO i_ao_1 = 1, n_ao
               DO i_ao_2 = 1, n_ao

                  ! Σ_λσ^R = sum_PR1νS1 M^G_λ0,νS1,PR1 * M^W_σR,νS1,PR1, where
                  ! M^G_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) G_µν^S2
                  ! M^W_σR,νS1,PR1 = sum_QR2 (σR νS1 | QR1-R2) W_PQ^R2 = M^W_σ0,νS1-R,PR1-R
                  sum_GW = SUM(M_G(i_ao_1, :, :)* &
                               M_W(i_cell_S1_minus_R, i_cell_R1_minus_R)%data_3(i_ao_2, :, :))
                  Sigma_R(i_ao_1, i_ao_2, i_cell_R) = Sigma_R(i_ao_1, i_ao_2, i_cell_R) &
                                                      + sign_Sigma*sum_GW
               END DO
            END DO

         END DO ! i_cell_R

         CALL timestop(handle2)

      END DO ! i_cell_R1

      CALL timestop(handle)

   END SUBROUTINE contract_to_Sigma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_t, i_task_local, ispin, &
 i_ao_1, i_ao_2, i_cell_R
      INTEGER, DIMENSION(3)                              :: cell_R
      INTEGER(KIND=int_8)                                :: n_double_in_3c
      LOGICAL                                            :: memory_fully_used
      REAL(KIND=dp)                                      :: t1, tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: G_occ_S, G_vir_S, Sigma_c_R_neg_tau, &
                                                            Sigma_c_R_pos_tau, W_R_tau
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_D, M_W
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c_G, store_int_3c_W

      CALL timeset(routineN, handle)

      t1 = m_walltime()
      n_double_in_3c = INT(0, KIND=int_8)

      CALL alloc_for_Sigma_c(M_D, M_W, store_int_3c_W, store_int_3c_G, G_occ_S, G_vir_S, &
                             Sigma_c_R_neg_tau, Sigma_c_R_pos_tau, W_R_tau, bs_env)

      memory_fully_used = .FALSE.

      ! Σ^c_λσ^R(iτ) = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) G^occ/vir_µν^S2(i|τ|) ]
      !                           [ sum_QR2 (σR νS1    | QR1-R2) Ŵ_PQ^R2(iτ)           ]
      DO i_t = 1, bs_env%num_time_freq_points
         DO ispin = 1, bs_env%n_spin

            t1 = m_walltime()

            tau = bs_env%imag_time_points(i_t)

            ! G^occ_µλ(i|τ|,k) = sum_n^occ C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k), τ < 0
            ! G^vir_µλ(i|τ|,k) = sum_n^vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k), τ > 0
            ! k-point k -> cell S: G^occ/vir_µλ^S(i|τ|) = sum_k w_k G^occ/vir_µλ(i|τ|,k) e^(ikS)
            CALL G_occ_vir(bs_env, tau, G_occ_S, ispin, occ=.TRUE., vir=.FALSE.)
            CALL G_occ_vir(bs_env, tau, G_vir_S, ispin, occ=.FALSE., vir=.TRUE.)

            ! get W_PQ^R(iτ)
            CALL fm_to_local_array(bs_env%fm_MWM_R_t(:, i_t), W_R_tau)

            Sigma_c_R_neg_tau(:,:,:) = 0.0_dp
            Sigma_c_R_pos_tau(:,:,:) = 0.0_dp

            ! JW 2del
            IF (bs_env%para_env%mepos == 0 ) THEN
               DO i_cell_R = 1, bs_env%nimages_scf_desymm
                  cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
                  DO i_ao_1 = 1, MIN(bs_env%n_ao, 30)
                     DO i_ao_2 = 1, MIN(bs_env%n_ao, 30)
                        PRINT *, 'tj =', i_t, 'λ =', i_ao_1, 'σ =', i_ao_2, 'R =', cell_R(1:2), &
                           'Gov_λσ^R(iτ_j mp 0) =', G_occ_S(i_ao_1, i_ao_2, i_cell_R), &
                              G_vir_S(i_ao_1, i_ao_2, i_cell_R)
                     END DO
                  END DO
               END DO
            END IF
           ! end 2del

            ! JW 2del
            IF (bs_env%para_env%mepos == 0 ) THEN
               DO i_cell_R = 1, bs_env%nimages_scf_desymm
                  cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
                  DO i_ao_1 = 1, MIN(bs_env%n_RI, 30)
                     DO i_ao_2 = 1, MIN(bs_env%n_RI, 30)
                        PRINT *, 'tj =', i_t, 'P =', i_ao_1, 'Q =', i_ao_2, 'R =', cell_R(1:2), &
                           'MWM_PQ^R(iτ_j mp 0) =', W_R_tau(i_ao_1, i_ao_2, i_cell_R)
                     END DO
                  END DO
               END DO
            END IF
           ! end 2del

            DO i_task_local = 1, bs_env%n_tasks_local

               ! M^W_σ0,νS1,PR1 = sum_QR2 ( σ0 νS1 | QR1-R2 ) W_QP^R2(iτ) for i_task_local
               CALL compute_3c_and_contract_W(M_W, W_R_tau, store_int_3c_W, i_task_local, &
                                              n_double_in_3c, memory_fully_used, bs_env, qs_env)

               ! M^Gocc_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) G^occ_µν^S2(i|τ|)
               ! Σ^c_λσ^R(iτ) = sum_PR1νS1 M^Gocc_λ0,νS1,PR1 * M^W_σR,νS1,PR1 for i_task_local, where
               ! ( τ < 0 )                                     M^W_σR,νS1,PR1 = M^W_σ0,νS1-R,PR1-R
               CALL contract_to_Sigma(Sigma_c_R_neg_tau, G_occ_S, M_W, store_int_3c_G, &
                                      i_task_local, n_double_in_3c, memory_fully_used, &
                                      bs_env, qs_env, occ=.TRUE., vir=.FALSE.)

               ! M^Gvir_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) G^vir_µν^S2(i|τ|)
               ! Σ^c_λσ^R(iτ) = sum_PR1νS1 M^Gocc_λ0,νS1,PR1 * M^W_σR,νS1,PR1 for i_task_local, where
               ! ( τ > 0 )                                     M^W_σR,νS1,PR1 = M^W_σ0,νS1-R,PR1-R
               CALL contract_to_Sigma(Sigma_c_R_pos_tau, G_vir_S, M_W, store_int_3c_G, &
                                      i_task_local, n_double_in_3c, memory_fully_used, &
                                      bs_env, qs_env, occ=.FALSE., vir=.TRUE.)

            END DO ! i_task_local: atom_P, atom_ν, i_cell_Delta_R

            CALL bs_env%para_env%sync()
            CALL bs_env%para_env%sum(Sigma_c_R_neg_tau)
            CALL bs_env%para_env%sum(Sigma_c_R_pos_tau)

            CALL local_array_to_fm(Sigma_c_R_neg_tau, bs_env%fm_Sigma_c_R_neg_tau(:, i_t))
            CALL local_array_to_fm(Sigma_c_R_pos_tau, bs_env%fm_Sigma_c_R_pos_tau(:, i_t))

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I10,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ,k=0) for time point ', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         ! JW 2del
         IF (bs_env%para_env%mepos == 0 .AND. i_t > 1) THEN
            DO i_cell_R = 1, bs_env%nimages_scf_desymm
               cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
               DO i_ao_1 = 1, MIN(bs_env%n_ao, 30)
                  DO i_ao_2 = 1, MIN(bs_env%n_ao, 30)
                     PRINT *, 'tj =', i_t, 'λ =', i_ao_1, 'σ =', i_ao_2, 'R =', cell_R(1:2), &
                        'Σ^c_λσ^R(iτ_j mp 0) =', Sigma_c_R_neg_tau(i_ao_1, i_ao_2, i_cell_R), &
                           Sigma_c_R_pos_tau(i_ao_1, i_ao_2, i_cell_R)
                  END DO
               END DO
            END DO
         END IF

         IF (bs_env%para_env%mepos == 0  .AND. i_t == 1) THEN
            DO i_cell_R = 1, bs_env%nimages_scf_desymm
               cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
               DO i_ao_1 = 1, bs_env%n_ao
                  DO i_ao_2 = 1, bs_env%n_ao
                     PRINT *, 'tj =', i_t, 'λ =', i_ao_1, 'σ =', i_ao_2, 'R =', cell_R(1:2), &
                        'Σ^c_λσ^R(iτ_j mp 0) =', Sigma_c_R_neg_tau(i_ao_1, i_ao_2, i_cell_R), &
                           Sigma_c_R_pos_tau(i_ao_1, i_ao_2, i_cell_R)
                  END DO
               END DO
            END DO
         END IF

         ! end 2del

         END DO ! ispin
      END DO ! i_t

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param M_D ...
!> \param M_W ...
!> \param store_int_3c_W ...
!> \param store_int_3c_G ...
!> \param G_occ_S ...
!> \param G_vir_S ...
!> \param Sigma_c_R_neg_tau ...
!> \param Sigma_c_R_pos_tau ...
!> \param W_R_tau ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE alloc_for_Sigma_c(M_D, M_W, store_int_3c_W, store_int_3c_G, G_occ_S, G_vir_S, &
                                Sigma_c_R_neg_tau, Sigma_c_R_pos_tau, W_R_tau, bs_env)

      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :)    :: M_D, M_W
      TYPE(data_3_type), ALLOCATABLE, DIMENSION(:, :, :) :: store_int_3c_W, store_int_3c_G
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: G_occ_S, G_vir_S, Sigma_c_R_neg_tau, &
                                                            Sigma_c_R_pos_tau, W_R_tau
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'alloc_for_Sigma_c'

      INTEGER                                            :: handle, n_ao, n_atom, n_RI, nimages_3c, &
                                                            nimages_scf

      CALL timeset(routineN, handle)

      nimages_3c = bs_env%nimages_3c
      nimages_scf = bs_env%nimages_scf_desymm
      n_atom = bs_env%n_atom
      n_ao = bs_env%n_ao
      n_RI = bs_env%n_RI

      ! storage for 3-index quantities
      ALLOCATE (M_D(nimages_3c, nimages_3c))
      ALLOCATE (M_W(nimages_3c, nimages_3c))
      ALLOCATE (store_int_3c_W(nimages_3c, nimages_3c, n_atom))
      ALLOCATE (store_int_3c_G(nimages_3c, nimages_3c, n_atom))

      ! G^occ_µλ(i|τ|)^S and G^vir_µλ(i|τ|)^S storage for single τ
      ALLOCATE (G_occ_S(n_ao, n_ao, nimages_scf))
      ALLOCATE (G_vir_S(n_ao, n_ao, nimages_scf))

      ! Σ^c_λσ^R(iτ) for τ < 0 (contraction with G^occ) and τ > 0 (contraction with G^virt)
      ! Sigma_c_R_neg_tau and Sigma_c_R_pos_tau is temporary storage for a single τ
      ALLOCATE (Sigma_c_R_neg_tau(n_ao, n_ao, nimages_scf))
      ALLOCATE (Sigma_c_R_pos_tau(n_ao, n_ao, nimages_scf))
      Sigma_c_R_neg_tau(:, :, :) = 0.0_dp
      Sigma_c_R_pos_tau(:, :, :) = 0.0_dp

      ! W_PQ^R(iτ) storage for single τ
      ALLOCATE (W_R_tau(n_RI, n_RI, nimages_scf))
      W_R_tau(:, :, :) = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE alloc_for_Sigma_c

   SUBROUTINE compute_QP_energies(bs_env, qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER                                            :: handle, ikp, ispin, j_t, i_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n, V_xc_ikp_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_freq, Sigma_c_ikp_n_time
      TYPE(cp_cfm_type)                                  :: cfm_ks_ikp, cfm_mo_coeff, &
                                                            cfm_Sigma_x_ikp, cfm_work_ikp

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mo_coeff, bs_env%fm_s_Gamma%matrix_struct)
      ALLOCATE (Sigma_x_ikp_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_ikp_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_ikp_n_freq(bs_env%n_ao, bs_env%num_time_freq_points, 2))

      DO ispin = 1, bs_env%n_spin

         DO ikp = 1, bs_env%nkp_bs

            ! 1. get C_µn(k)
            CALL cp_fm_to_cfm(bs_env%fm_mo_coeff_kp(ikp, ispin, 1), &
                              bs_env%fm_mo_coeff_kp(ikp, ispin, 2), cfm_mo_coeff)

            ! 2. Σ^x_µν(k) = sum_R Σ^x_µν^R e^ikR
            !    Σ^x_nn(k) = sum_µν C^*_µn(k) Σ^x_µν(k) C_νn(k) 
            CALL trafo_to_k_and_nn(bs_env%fm_Sigma_x_R, Sigma_x_ikp_n, cfm_mo_coeff, bs_env, ikp)

            ! 3. Σ^c_µν(k,+/-i|τ_j|) = sum_R Σ^c_µν^R(+/-i|τ_j|) e^ikR
            !    Σ^c_nn(k,+/-i|τ_j|) = sum_µν C^*_µn(k) Σ^c_µν(k,+/-i|τ_j|) C_νn(k)
            DO j_t = 1, bs_env%num_time_freq_points
                CALL trafo_to_k_and_nn(bs_env%fm_Sigma_c_R_pos_tau(:, j_t), &
                                       Sigma_c_ikp_n_time(:, j_t, 1), cfm_mo_coeff, bs_env, ikp)
                CALL trafo_to_k_and_nn(bs_env%fm_Sigma_c_R_neg_tau(:, j_t), &  
                                       Sigma_c_ikp_n_time(:, j_t, 2), cfm_mo_coeff, bs_env, ikp)
            END DO

            ! JW 2del
            IF(ikp == 1 .AND. bs_env%para_env%mepos == 0) THEN
              DO j_t = 1, bs_env%num_time_freq_points
                PRINT *, 'j_t =', j_t, 'Σ_11(+/-τ_j,k_1) =', Sigma_c_ikp_n_time(1, j_t, 1:2)
              END DO
              DO j_t = 1, bs_env%num_time_freq_points
                IF(bs_env%n_occ(1) > 1) THEN
                PRINT *, 'j_t =', j_t, 'Σ_22(+/-τ_j,k_1) =', Sigma_c_ikp_n_time(2, j_t, 1:2)
                END IF
              END DO
              DO j_t = 1, bs_env%num_time_freq_points
                IF(bs_env%n_occ(1) > 2) THEN
                PRINT *, 'j_t =', j_t, 'Σ_33(+/-τ_j,k_1) =', Sigma_c_ikp_n_time(3, j_t, 1:2)
                END IF
              END DO
            END IF

            ! 4. Σ^c_nn(k_i,iω) = ∫ from -∞ to ∞ dτ e^-iωτ Σ^c_nn(k_i,iτ)
            CALL time_to_freq(bs_env, Sigma_c_ikp_n_time, Sigma_c_ikp_n_freq, ispin)

            ! 5. Analytic continuation Σ^c_nn(k_i,iω) -> Σ^c_nn(k_i,ϵ) and
            !    ϵ_nk_i^GW = ϵ_nk_i^DFT + Σ^c_nn(k_i,ϵ) + Σ^x_nn(k_i) - v^xc_nn(k_i)
            CALL analyt_conti_and_print(bs_env, Sigma_c_ikp_n_freq, Sigma_x_ikp_n, &
                                        bs_env%v_xc_n(:, ikp, ispin), &
                                        bs_env%eigenval_scf(:, ikp, ispin), ikp, &
                                        bs_env%nkp_bs, ispin)

         END DO ! ikp

      END DO ! ispin

      CALL get_VBM_CBM_bandgaps(bs_env)

      CALL cp_cfm_release(cfm_mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies


   SUBROUTINE trafo_to_k_and_nn(fm_rs, array_ikp_n, cfm_mo_coeff, bs_env, ikp)
      TYPE(cp_fm_type), DIMENSION(:)    :: fm_rs
      REAL(KIND=dp), DIMENSION(:)           :: array_ikp_n
      TYPE(cp_cfm_type)                                  :: cfm_mo_coeff
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'trafo_to_k_and_nn'

      INTEGER                                            :: handle, ikp, ispin, j_t, n_ao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n, V_xc_ikp_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_freq, Sigma_c_ikp_n_time
      TYPE(cp_cfm_type)                                  :: cfm_ks_ikp, cfm_mo_coeff_ikp, cfm_tmp, &
                                                            cfm_Sigma_x_ikp, cfm_work_ikp, cfm_ikp
      TYPE(cp_fm_type)                                  :: fm_ikp_re, fm_ikp_im
      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_ikp, cfm_mo_coeff%matrix_struct)       
      CALL cp_cfm_create(cfm_tmp, cfm_mo_coeff%matrix_struct)
      CALL cp_fm_create(fm_ikp_re, cfm_mo_coeff%matrix_struct)

      ! Σ_µν(k_i) = sum_R e^ik_iR Σ_µν^R
      CALL fm_trafo_rs_to_ikp(cfm_ikp, fm_rs, bs_env%kpoints_DOS, ikp)

      n_ao = bs_env%n_ao

      ! Σ_nm(k_i) = sum_µν C^*_µn(k_i) Σ_µν(k_i) C_νn(k_i)
      CALL parallel_gemm('N', 'N', n_ao, n_ao, n_ao, z_one, cfm_ikp, cfm_mo_coeff, z_zero, cfm_tmp)
      CALL parallel_gemm('C', 'N', n_ao, n_ao, n_ao, z_one, cfm_mo_coeff, cfm_tmp, z_zero, cfm_ikp)

      ! get Σ_nn(k_i) which is a real quantity as Σ^x and Σ^c(iτ) is Hermitian
      CALL cp_cfm_to_fm(cfm_ikp, fm_ikp_re)
      CALL cp_fm_get_diag(fm_ikp_re, array_ikp_n)

      CALL cp_cfm_release(cfm_ikp)
      CALL cp_cfm_release(cfm_tmp)
      CALL cp_fm_release(fm_ikp_re)

      CALL timestop(handle)

   END SUBROUTINE trafo_to_k_and_nn

END MODULE gw_small_cell_full_kp
