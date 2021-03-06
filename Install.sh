# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

  cp types.h ..

  cp bond_hf_harmonic.h ..
  cp bond_hf_harmonic.cpp ..

  cp atom_vec_ch_bond.cpp .. 
  cp atom_vec_ch_bond.h .. 

  cp create_hf_atoms.cpp ..
  cp create_hf_atoms.h ..

  cp fix_injection_update.cpp ..
  cp fix_injection_update.h ..

  cp fix_lubrication_update_BC.cpp ..
  cp fix_lubrication_update_BC.h ..

  cp seedcrack1.cpp ..
  cp seedcrack1.h ..


  cp bond_rock_rock.cpp ..
  cp bond_rock_rock.h ..

  cp bond_channel_channel.cpp ..
  cp bond_channel_channel.h ..

  cp bond_rock_channel.cpp ..
  cp bond_rock_channel.h ..

  cp bond_channel_rock.cpp ..
  cp bond_channel_rock.h ..

  cp deformation.cpp ..
  cp deformation.h ..

  cp dump_hf_image.cpp ..
  cp dump_hf_image.h ..

  cp dump_hf_image2.cpp ..
  cp dump_hf_image2.h ..

  cp check_list.cpp ..
  cp check_list.h ..

  cp compute_channel_pressure.cpp ..
  cp compute_channel_pressure.h ..

  cp dump_ch_local.cpp ..
  cp dump_ch_local.h ..

  cp delta_cutoff_randomness2.cpp .. 
  cp delta_cutoff_randomness2.h .. 

  cp delta_cutoff_randomness3.cpp .. 
  cp delta_cutoff_randomness3.h .. 

elif (test $1 = 0) then

  rm ../types.h

  rm ../bond_hf_harmonic.h
  rm ../bond_hf_harmonic.cpp

  rm ../atom_vec_ch_bond.cpp 
  rm ../atom_vec_ch_bond.h

  rm ../create_hf_atoms.cpp
  rm ../create_hf_atoms.h

  rm ../fix_injection_update.cpp
  rm ../fix_injection_update.h
 
  rm ../fix_lubrication_update_BC.cpp
  rm ../fix_lubrication_update_BC.h

  rm ../seedcrack1.cpp
  rm ../seedcrack1.h

  rm ../bond_rock_rock.cpp
  rm ../bond_rock_rock.h

  rm ../bond_channel_channel.cpp
  rm ../bond_channel_channel.h

  rm ../bond_rock_channel.cpp
  rm ../bond_rock_channel.h

  rm ../bond_channel_rock.cpp
  rm ../bond_channel_rock.h

  rm ../deformation.cpp
  rm ../deformation.h

  rm ../dump_hf_image.cpp
  rm ../dump_hf_image.h

  rm ../dump_hf_image2.cpp
  rm ../dump_hf_image2.h

  rm ../check_list.cpp
  rm ../check_list.h

  rm ../compute_channel_pressure.cpp
  rm ../compute_channel_pressure.h

  rm ../dump_ch_local.cpp
  rm ../dump_ch_local.h

  rm ../delta_cutoff_randomness2.cpp 
  rm ../delta_cutoff_randomness2.h 

  rm ../delta_cutoff_randomness3.cpp 
  rm ../delta_cutoff_randomness3.h 


fi
