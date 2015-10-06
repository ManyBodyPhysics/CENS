MODULE chp_presets

  use constants
  use idaho_chiral_potential
  
CONTAINS
    
  subroutine chp_preset(preset)
  
    character(len=100), INTENT(IN) :: preset
    character(len=42) :: title
    integer :: print_to_screen=6
    
    if (type_of_pot /= 'chiral-potential') RETURN
    
    SELECT CASE(preset)
    CASE('nnloopt')
       title = 'nnloopt. PRL 110(19), 192502 (2013)'
       call chp_preset_nnloopt
    CASE('nnlo400')
       title = 'nnlo400:LAM=400/SFR=700'
       call chp_preset_nnlo400
    CASE('nnlo450')
       title = 'nnlo450:LAM=450/SFR=700'
       call chp_preset_nnlo450
    CASE('nnlo500')
       title = 'nnlo500:LAM=500/SFR=700'
       call chp_preset_nnlo500
    CASE('nnlo550')
       title = 'nnlo550:LAM=550/SFR=700'
       call chp_preset_nnlo550
    CASE('nnlo600')
       title = 'nnlo600:LAM=600/SFR=700'
       call chp_preset_nnlo600
    CASE('nnlo650')
       title = 'nnlo650:LAM=650/SFR=700'
       call chp_preset_nnlo650
    CASE('nnlo700')
       title = 'nnlo700:LAM=700/SFR=700'
       call chp_preset_nnlo700
    CASE('nnlo750')
       title = 'nnlo750:LAM=750/SFR=700'
       call chp_preset_nnlo750
    CASE('idaho-n3lo500')
       title = 'idaho-n3lo500: Rup''s N3LO LAM=500/DR'
       call chp_preset_idaho_n3lo500
    CASE DEFAULT
       write(*,*) 'chp_preset: unknown potential'
    END SELECT
    
    write(*,*) 'using chiral potential: ', title
    write(*,*) 'from chiral-twobody-potential.f90 '
    write(*,*) 'see definition in output_run file'
    
    call chp_get_mass_nucleon(p_mass)
    p_massave = p_mass(0)
    
    call chp_print_constants(title,print_to_screen)
    
  end subroutine chp_preset
  
  ! --- BEGIN AUTO-GENERATED --- 140411
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/NNLOopt.ini +potential_name=nnloopt
  subroutine chp_preset_nnloopt
    use idaho_chiral_potential
    
    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(500.000D0)

    call chp_set_c1(  -0.918639528734720D0)
    call chp_set_c3(  -3.888687492763241D0)
    call chp_set_c4(   4.310327160829740D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.151366037203108D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.158434176622812D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.152141088236679D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.158434176622812D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.151764745900691D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.158434176622812D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.404021944134705D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.263390763475578D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.417045542055649D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.782658499975205D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.928384662662304D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.618141419047458D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.677808511406356D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_2PE_CSB_correct_mass(1)

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnloopt
  ! --- END AUTO-GENERATED ---
  
  ! --- BEGIN AUTO-GENERATED --- 140410
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut400sfr700-potential_125.ini +potential_name=nnlo400
  subroutine chp_preset_nnlo400
    use idaho_chiral_potential
    
    implicit none
    
    call initialize_chiral_potential
    
    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))
    
    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(400.000D0)

    call chp_set_c1(  -0.911687940175766D0)
    call chp_set_c3(  -3.880267525993659D0)
    call chp_set_c4(   5.100215517413014D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.152629075685936D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.182297517383481D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.153468020417515D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.182297517383481D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.153074751601327D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.182297517383481D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.455478237837730D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.226304720793696D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.610215953642314D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.897883207536673D0) ! C_3P1
    call chp_set_NLO_contact(5,    1.068865994908072D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.747836527329264D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.678070210121440D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants
    
  end subroutine chp_preset_nnlo400
! --- END AUTO-GENERATED ---

! --- BEGIN AUTO-GENERATED ---
!   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut450sfr700-potential_125.ini +potential_name=nnlo450

subroutine chp_preset_nnlo450
    use idaho_chiral_potential

    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(450.000D0)

    call chp_set_c1(  -0.910294817329351D0)
    call chp_set_c3(  -3.880687659323964D0)
    call chp_set_c4(   4.670920622150460D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.152035456748769D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.169539566807816D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.152827403337367D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.169539566807816D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.152472577107566D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.169539566807816D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.431098288036233D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.215167443048370D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.466918210213105D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.850349848943376D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.987574362681558D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.681421332367734D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.673182679211219D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo450
  ! --- END AUTO-GENERATED ---
  ! --- BEGIN AUTO-GENERATED ---
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut500sfr700-potential_125.ini +potential_name=nnlo500
  
  subroutine chp_preset_nnlo500
    use idaho_chiral_potential
    
    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(500.000D0)

    call chp_set_c1(  -0.919407456831721D0)
    call chp_set_c3(  -3.889838481170248D0)
    call chp_set_c4(   4.307367469348735D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.151363639522344D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.158481250438791D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.152152631085613D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.158481250438791D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.151804824990968D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.158481250438791D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.404312353608718D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.265789779817624D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.414829078808030D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.779984839303950D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.927937115539244D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.618550395799233D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.673470419600444D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_2PE_CSB_correct_mass(1)

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo500
  ! --- END AUTO-GENERATED ---
  ! --- BEGIN AUTO-GENERATED ---
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut550sfr700-potential_125.ini +potential_name=nnlo550
  
  subroutine chp_preset_nnlo550
    use idaho_chiral_potential

    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(550.000D0)

    call chp_set_c1(  -0.906302676626195D0)
    call chp_set_c3(  -3.897385334719676D0)
    call chp_set_c4(   3.906282430405821D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.150672782581269D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.146778629499074D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.151623708436515D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.146778629499074D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.151215788399478D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.146778629499074D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.389653893342561D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.325329841680243D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.386120512442197D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.684247439232355D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.838995783502866D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.562661202867282D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.674440897464630D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo550
  ! --- END AUTO-GENERATED ---
  ! --- BEGIN AUTO-GENERATED ---
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut600sfr700-potential_125.ini +potential_name=nnlo600
  
  subroutine chp_preset_nnlo600
    use idaho_chiral_potential

    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(600.000D0)

    call chp_set_c1(  -0.904098177403228D0)
    call chp_set_c3(  -3.900923127946868D0)
    call chp_set_c4(   3.591167936015824D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.149108165110814D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.137586184728017D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.150166048702156D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.137586184728017D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.149707974692694D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.137586184728017D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.362005506385191D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.511206053190680D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.383185320481828D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.522273640822317D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.778278220846694D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.517472890926010D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.682867912136692D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo600
  ! --- END AUTO-GENERATED ---
  ! --- BEGIN AUTO-GENERATED ---
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut650sfr700-potential_125.ini +potential_name=nnlo650
  
  subroutine chp_preset_nnlo650
    use idaho_chiral_potential

    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(650.000D0)

    call chp_set_c1(  -0.927701587626971D0)
    call chp_set_c3(  -3.886122146126503D0)
    call chp_set_c4(   3.395513351863698D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.147062207891242D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.129226506464466D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.148345486147834D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.129226506464466D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.147787847438749D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.129226506464466D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.345181161539688D0) ! C_1S0
    call chp_set_NLO_contact(2,    2.104925093661719D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.400381214737342D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.195844901116646D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.725225140031640D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.485207168344209D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.680445992619058D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo650
  ! --- END AUTO-GENERATED ---
  ! --- BEGIN AUTO-GENERATED ---
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut700sfr700-potential_125.ini +potential_name=nnlo700
  
  subroutine chp_preset_nnlo700
    use idaho_chiral_potential
    
    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(700.000D0)

    call chp_set_c1(  -0.914211988510583D0)
    call chp_set_c3(  -3.888884071668157D0)
    call chp_set_c4(   3.285545596366145D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.143841733895033D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.122136548151004D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.145439634819476D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.122136548151004D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.144732495927184D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.122136548151004D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.338000463643345D0) ! C_1S0
    call chp_set_NLO_contact(2,    5.966123165645889D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.435388650413494D0) ! C_1P1
    call chp_set_NLO_contact(4,    0.819408982995883D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.696102531702483D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.464829987199691D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.672824002909096D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo700
  ! --- END AUTO-GENERATED ---
  ! --- BEGIN AUTO-GENERATED ---
  !   CMD: ./make_chp_init_code ./local_input/potential_definitions/nnlo-input/n2lo_cut750sfr700-potential_125.ini +potential_name=nnlo750
  
  subroutine chp_preset_nnlo750
    use idaho_chiral_potential
    
    implicit none

    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9184000000D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(NNLO)
    call chp_set_reg("SF", 700.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(750.000D0)

    call chp_set_c1(  -0.911597474931629D0)
    call chp_set_c3(  -3.897878976933375D0)
    call chp_set_c4(   3.181806412940659D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.139386737492394D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.116059224500765D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.141428604836242D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.116059224500765D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.140501064845088D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.116059224500765D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.331677562259216D0) ! C_1S0
    call chp_set_NLO_contact(2,   28.874052342211318D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.582591203161031D0) ! C_1P1
    call chp_set_NLO_contact(4,    1.927270998797956D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.663274299417975D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.447246360982541D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.667996190893353D0) ! C_3P2

    call chp_set_1PE_reg_par(3.0D0)
    call chp_set_2PE_reg_par(3.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 3.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 3.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 3.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 3.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 3.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 3.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 3.0D0) ! C_3P2

    call chp_set_units_and_derive_constants

  end subroutine chp_preset_nnlo750
  ! --- END AUTO-GENERATED ---
    
! --- BEGIN AUTO-GENERATED ---
!   CMD: ./make_chp_init_code ./local_input/potential_definitions/idaho-n3lo-potential.ini +potential_name=idaho_n3lo500

  subroutine chp_preset_idaho_n3lo500
    use idaho_chiral_potential
    
    implicit none
    
    call initialize_chiral_potential

    ! (proton, neutron, nucleon)
    call chp_set_mass_nucleon((/938.2720000000D0, 939.5653000000D0, 938.9182046406D0/))
    ! (pi-, pi, pi+)
    call chp_set_mass_pion((/139.5702000000D0, 134.9766000000D0, 139.5702000000D0/))

    call chp_set_chiral_order(N3LO)
    call chp_set_reg("DR", 0.000D0)
    call chp_set_itope("EM")
    call chp_set_contact_format("PW")

    call chp_set_gA(1.2900D0)
    call chp_set_fpi(92.4000D0)
    call chp_set_fine_structure(0.007297352570000D0)

    call chp_set_Lambda(500.000D0)

    call chp_set_c1(  -0.810000000000000D0)
    call chp_set_c3(  -3.200000000000000D0)
    call chp_set_c4(   5.400000000000000D0)
    call chp_set_c2           (   2.800000000000000D0)
    call chp_set_d1_plus_d2   (   3.060000000000000D0)
    call chp_set_d3           (  -3.270000000000000D0)
    call chp_set_d5           (   0.450000000000000D0)
    call chp_set_d14_minus_d15(  -5.650000000000000D0)

    call chp_set_CIB_LO_contact(1, -1,   -0.145286000000000D0) ! Ct_1S0pp
    call chp_set_CIB_LO_contact(2, -1,   -0.118972496000000D0) ! Ct_3S1pp
    call chp_set_CIB_LO_contact(1,  0,   -0.147167000000000D0) ! Ct_1S0np
    call chp_set_CIB_LO_contact(2,  0,   -0.118972496000000D0) ! Ct_3S1np
    call chp_set_CIB_LO_contact(1,  1,   -0.146285000000000D0) ! Ct_1S0nn
    call chp_set_CIB_LO_contact(2,  1,   -0.118972496000000D0) ! Ct_3S1nn

    call chp_set_NLO_contact(1,    2.380000000000000D0) ! C_1S0
    call chp_set_NLO_contact(2,    1.487000000000000D0) ! C_3P0
    call chp_set_NLO_contact(3,    0.656000000000000D0) ! C_1P1
    call chp_set_NLO_contact(4,   -0.630000000000000D0) ! C_3P1
    call chp_set_NLO_contact(5,    0.760000000000000D0) ! C_3S1
    call chp_set_NLO_contact(6,    0.826000000000000D0) ! C_3S1-3D1
    call chp_set_NLO_contact(7,   -0.538000000000000D0) ! C_3P2

    call chp_set_N3LO_contact( 1,   -2.545000000000000D0) ! Dh_1S0
    call chp_set_N3LO_contact( 2,  -16.000000000000000D0) ! D_1S0
    call chp_set_N3LO_contact( 3,    0.245000000000000D0) ! D_3P0
    call chp_set_N3LO_contact( 4,    5.250000000000000D0) ! D_1P1
    call chp_set_N3LO_contact( 5,    2.350000000000000D0) ! D_3P1
    call chp_set_N3LO_contact( 6,    7.000000000000001D0) ! Dh_3S1
    call chp_set_N3LO_contact( 7,    6.550000000000000D0) ! D_3S1
    call chp_set_N3LO_contact( 8,   -2.800000000000000D0) ! D_3D1
    call chp_set_N3LO_contact( 9,    2.250000000000000D0) ! Dh_3S1-3D1
    call chp_set_N3LO_contact(10,    6.610000000000000D0) ! D_3S1-3D1
    call chp_set_N3LO_contact(11,   -1.770000000000000D0) ! D_1D2
    call chp_set_N3LO_contact(12,   -1.460000000000000D0) ! D_3D2
    call chp_set_N3LO_contact(13,    2.295000000000000D0) ! D_3P2
    call chp_set_N3LO_contact(14,   -0.465000000000000D0) ! D_3P2-3F2
    call chp_set_N3LO_contact(15,    5.660000000000000D0) ! D_3D3

    call chp_set_1PE_reg_par(4.0D0)
    call chp_set_2PE_reg_par(2.0D0)
    call chp_set_LO_contact_reg_par(1, 3.0D0) ! Ct_1S0
    call chp_set_LO_contact_reg_par(2, 3.0D0) ! Ct_3S1
    call chp_set_NLO_contact_reg_par(1, 2.0D0) ! C_1S0
    call chp_set_NLO_contact_reg_par(2, 2.0D0) ! C_3P0
    call chp_set_NLO_contact_reg_par(3, 2.0D0) ! C_1P1
    call chp_set_NLO_contact_reg_par(4, 2.0D0) ! C_3P1
    call chp_set_NLO_contact_reg_par(5, 2.0D0) ! C_3S1
    call chp_set_NLO_contact_reg_par(6, 2.0D0) ! C_3S1-3D1
    call chp_set_NLO_contact_reg_par(7, 2.0D0) ! C_3P2
    call chp_set_N3LO_contact_reg_par(1, 2.0D0) ! Dh_1S0
    call chp_set_N3LO_contact_reg_par(2, 2.0D0) ! D_1S0
    call chp_set_N3LO_contact_reg_par(3, 3.0D0) ! D_3P0
    call chp_set_N3LO_contact_reg_par(4, 2.0D0) ! D_1P1
    call chp_set_N3LO_contact_reg_par(5, 4.0D0) ! D_3P1
    call chp_set_N3LO_contact_reg_par(6, 2.0D0) ! Dh_3S1
    call chp_set_N3LO_contact_reg_par(7, 2.0D0) ! D_3S1
    call chp_set_N3LO_contact_reg_par(8, 2.0D0) ! D_3D1
    call chp_set_N3LO_contact_reg_par(9, 2.0D0) ! Dh_3S1-3D1
    call chp_set_N3LO_contact_reg_par(10, 2.0D0) ! D_3S1-3D1
    call chp_set_N3LO_contact_reg_par(11, 4.0D0) ! D_1D2
    call chp_set_N3LO_contact_reg_par(12, 2.0D0) ! D_3D2
    call chp_set_N3LO_contact_reg_par(13, 2.0D0) ! D_3P2
    call chp_set_N3LO_contact_reg_par(14, 4.0D0) ! D_3P2-3F2
    call chp_set_N3LO_contact_reg_par(15, -1.0D0) ! D_3D3

    call chp_set_units_and_derive_constants
    
  end subroutine chp_preset_idaho_n3lo500
! --- END AUTO-GENERATED ---


END MODULE chp_presets
