/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(injection/update,FixInjectionUpdate)

#else

#ifndef LMP_FIX_INJECTION_UPDATE_H
#define LMP_FIX_INJECTION_UPDATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixInjectionUpdate : public Fix {
 public:
  FixInjectionUpdate(class LAMMPS *, int, char **);
  ~FixInjectionUpdate();
  //  void init();
  //  void init_list(int, class NeighList *);
  int setmask();
  void injection_process();

  void channel_update();
  int find_channel_atom(int, int);

  void pre_force(int);
  //void final_integrate();
 private:

  double injection_x, injection_y,injection_z;
  double injection_rate;
};

}

#endif
#endif
