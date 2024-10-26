### Mirror Repo: dyGiLa -- dynamical simulation of GL Effective Theory
HPC software for simulating the order parameter dynamics of superfluid Helium 3 with Time-Dependent Ginzburg-Landau effective theory.

dyGiLa started as a descendant of `onsim` developed by Dr. Asier Lopez-Eiguren. 
Same as `onsim`, dyGiLa uses `HILA` as its CFT simulation framework (Kari Rummukainen et al): https://cft-hy.github.io/HILA.home
After couple of years' develpment, Dr. Kuang Zhang has added parallel simulaton data stream feature and strong coupling corrections functions (material parameters) into dyGiLa.
Moreover, implemented by Dr. Asier Lopez-Eiguren, and generalized by Dr. Kuang Zhang, GL EOM with Langevin noise i.e., GL-Langevin equation has been introduced into dyGiLa with many different temperature profiles. 

### dyGiLa source tree
~~~ shellscript
 |-pario
 | |-src
 | | |-xml.cpp
 | | |-xdmf.cpp
 | | |-init.cpp
 | | |-shutdown.cpp
 | | |-utilities
 | | | |-actions_massCurrent.cpp
 | | | |-mesh_AMatrix.cpp
 | | | |-actions_spinCurrent.cpp
 | | | |-mesh_massCurrent.cpp
 | | | |-actions_printTree.cpp
 | | | |-mesh_addGhost_verify.cpp
 | | | |-mesh_gapA_FEDensity.cpp
 | | | |-mesh_spinCurrent.cpp
 | | | |-actions_AMatrix.cpp
 | | | |-actions_gapA_FEDensity.cpp
 | | | |-mesh.cpp
 | | |-pstream.cpp
 | |-inc
 | | |-pario.hpp
 | |-pario_conf.mk
 |-paras_conf
 | |-parameters_computed.txt
 | |-parameters_fixed.txt
 | |-sim_T0.txt
 | |-sim_params.txt
 | |-parameters_interpolated.txt
 | |-sim_config_dyGiLa-Langevin-quench.txt
 | |-sim_config_pario_Temperature_field.txt
 |-main.cpp
 |-Makefile
 |-README.md
 |-.gitignore
 |-#README.md#
 |-matep
 | |-matep_conf.mk
 | |-src
 | | |-matep.cpp
 | | |-matep_utils.cpp
 | |-inc
 | | |-matep.hpp
 |-glsol
 | |-src
 | | |-initialize
 | | | |-glsol_initialize.cpp
 | | | |-glsol_initialize_p.cpp
 | | | |-glsol_initialize_T.cpp
 | | |-allocate.cpp
 | | |-utilities
 | | | |-write_energies.cpp
 | | | |-write_positions.cpp
 | | | |-write_phases.cpp
 | | | |-point_params.cpp
 | | | |-hot_bloob.cpp
 | | | |-write_moduli.cpp
 | | |-next
 | | | |-next_bath_UniT_quench.cpp
 | | | |-next.cpp
 | | | |-next_T.cpp
 | | | |-next_bath.cpp
 | | | |-next_bath_UniT_quench.cpp~
 | |-inc
 | | |-dyGiLa_config.hpp
 | | |-glsol.hpp
 | |-glsol_conf.mk
~~~
