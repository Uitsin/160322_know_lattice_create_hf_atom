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

/* ----------------------------------------------------------------------
   this code is used to find the neighbors of channel elements (square lattice)
------------------------------------------------------------------------- */

#include "string.h"
#include "math.h"
#include "deformation.h"
#include "group.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "stdlib.h" 
#include "math.h" 
#include "mpi.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "memory.h"
#include "domain.h"
#include "timer.h"
#define TINY  1.e-3 ;
#define FAKE_INT_VALUE  -991 ;

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

Deformation::Deformation(LAMMPS *lmp) : Pointers(lmp) {}

void Deformation::command(int narg, char **arg)
{

  if (narg < 3) error->all(FLERR,"Illegal deformation command"); 
 

  double  xfactor = atof(arg[0]);
  double  yfactor = atof(arg[1]);
  double  zfactor = atof(arg[2]);

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int n ;
  
  int **bond_atom = atom -> bond_atom;
  int *tag = atom->tag;
  /*
  for (n = 0; n < nlocal; n++){
    if (tag[n] ==681){ fprintf(screen, "b==bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]); }
  }
  */
  for (n = 0; n < nlocal; n++){
    x[n][0] = xfactor* x[n][0];
    x[n][1] = yfactor* x[n][1];
    x[n][2] = zfactor* x[n][2];
    //    fprintf(screen,"the %d-th atom. xfactor yfactor %f %f x y %f %f \n", n,xfactor,yfactor, x[n][0],x[n][1]);
  }

  /*

  for (n = 0; n < nlocal; n++){
    if (tag[n] ==681){
      fprintf(screen, "d==bond_rock1,rock2, rock3, rock4 %d %d %d %d\n",bond_atom[n][0],bond_atom[n][1],bond_atom[n][2],bond_atom[n][3]);
    }
  }
  */
  timer->stamp();
  comm->forward_comm();
  timer->stamp(TIME_COMM);
}

