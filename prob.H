#ifndef PROB_H
#define PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "Constants.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  const amrex::Real xlo,
  const amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector,
  const ProbParmDevice& prob_parm)
{
  if (prob_parm.pmf_do_average) {
    int lo_loside = 0;
    int lo_hiside = 0;
    int hi_loside = 0;
    int hi_hiside = 0;
    if (xlo < prob_parm.d_pmf_X[0]) {
      lo_loside = 0;
      lo_hiside = 0;
    }
    if (xhi < prob_parm.d_pmf_X[0]) {
      hi_loside = 0;
      hi_hiside = 0;
    }
    if (xlo > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      lo_loside = prob_parm.pmf_N - 1;
      lo_hiside = prob_parm.pmf_N - 1;
    }
    if (xhi > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      hi_loside = prob_parm.pmf_N - 1;
      hi_hiside = prob_parm.pmf_N - 1;
    }
    if (lo_loside == 0) {
      for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if ((xlo > prob_parm.d_pmf_X[i]) && (xlo < prob_parm.d_pmf_X[i + 1])) {
          lo_loside = i;
          lo_hiside = i + 1;
        }
      }
    }
    if (hi_loside == 0) {
      for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if ((xhi > prob_parm.d_pmf_X[i]) && (xhi < prob_parm.d_pmf_X[i + 1])) {
          hi_loside = i;
          hi_hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < prob_parm.pmf_M; j++) {
      amrex::Real x1 = prob_parm.d_pmf_X[lo_loside];
      amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + lo_loside];
      amrex::Real x2 = prob_parm.d_pmf_X[lo_hiside];
      amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + lo_hiside];
      amrex::Real dydx = 0.0;
      if (lo_loside == lo_hiside) {
        dydx = 0.0;
      } else {
        dydx = (y2 - y1) / (x2 - x1);
      }
      amrex::Real ylo = y1 + dydx * (xlo - x1);
      amrex::Real yhi = 0.0;
      if (lo_loside == hi_loside) {
        yhi = y1 + dydx * (xhi - x1);
        y_vector[j] = 0.5 * (ylo + yhi);
      } else {
        amrex::Real sum = (x2 - xlo) * 0.5 * (ylo + y2);
        x1 = prob_parm.d_pmf_X[hi_loside];
        y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hi_loside];
        x2 = prob_parm.d_pmf_X[hi_hiside];
        y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hi_hiside];
        if (hi_loside == hi_hiside) {
          dydx = 0.0;
        } else {
          dydx = (y2 - y1) / (x2 - x1);
        }
        yhi = y1 + dydx * (xhi - x1);
        sum += (xhi - x1) * 0.5 * (yhi + y1);
        for (int k = lo_hiside; k < hi_loside - 1; k++) {
          sum += (prob_parm.d_pmf_X[k + 1] - prob_parm.d_pmf_X[k]) * 0.5 *
                 (prob_parm.d_pmf_Y[prob_parm.pmf_N * j + k] +
                  prob_parm.d_pmf_Y[prob_parm.pmf_N * j + k + 1]);
        }
        y_vector[j] = sum / (xhi - xlo);
      }
    }
  } else {
    amrex::Real xmid = 0.5 * (xlo + xhi);
    int loside = -1;
    int hiside = -1;
    if (xmid < prob_parm.d_pmf_X[0]) {
      loside = 0;
      hiside = 0;
    }
    if (xmid > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      loside = prob_parm.pmf_N - 1;
      hiside = prob_parm.pmf_N - 1;
    }
    if (loside == -1) {
      for (int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if (
          (xmid >= prob_parm.d_pmf_X[i]) &&
          (xmid <= prob_parm.d_pmf_X[i + 1])) {
          loside = i;
          hiside = i + 1;
        }
      }
    }
    for (int j = 0; j < prob_parm.pmf_M; j++) {
      const amrex::Real x1 = prob_parm.d_pmf_X[loside];
      const amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + loside];
      const amrex::Real x2 = prob_parm.d_pmf_X[hiside];
      const amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hiside];
      amrex::Real dydx = 0.0;
      if (loside == hiside) {
        dydx = 0.0;
      } else {
        dydx = (y2 - y1) / (x2 - x1);
      }
      y_vector[j] = y1 + dydx * (xlo - x1);
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
   // Geometry
  
   const amrex::Real* prob_lo = geomdata.ProbLo();
   const amrex::Real* prob_hi = geomdata.ProbHi();
   const amrex::Real* dx = geomdata.CellSize();

   const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
   const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
   const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

   const amrex::Real L_y = prob_hi[1] - prob_lo[1];

  //amrex::Real x_non_dim = x/prob_parm.delta_half;

  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  for (int n = 0; n < NUM_SPECIES; n++)
    massfrac[n] = prob_parm.massfrac[n];
    
  //state(i, j, k, USOOT) = prob_parm.P_amb;
  
  //Random density perturbations
  if (x<= prob_parm.x_ZND - prob_parm.A * std::sin(2*prob_parm.n*M_PI*y/L_y)) // || (x-0.25)*(x-0.25)+(y-0.25)*(y-0.25)-0.05*0.05<=0 || (x-0.25)*(x-0.25)+(y-0.5)*(y-0.5)-0.05*0.05<=0 || (x-0.25)*(x-0.25)+(y-0.75)*(y-0.75)-0.05*0.05<=0 ) //*std::sin(prob_parm.n*M_PI*z)) //3d sinusoidal
  {

     //assign
     state(i, j, k, URHO) = prob_parm.rho_amb;
     state(i, j, k, UMX) =  1.0*prob_parm.rho_amb * prob_parm.D_CJ;
     state(i, j, k, UMY) = 0.0;
     state(i, j, k, UMZ) = 0.0;
     state(i, j, k, UEINT) =  prob_parm.rho_amb *  prob_parm.eint_amb;
     state(i, j, k, UEDEN) =  prob_parm.rho_amb * ( prob_parm.eint_amb +  1.0 * prob_parm.D_CJ * prob_parm.D_CJ );
     state(i, j, k, UTEMP) = prob_parm.T_amb;
     for (int n = 0; n < NUM_SPECIES; n++)
        state(i, j, k, UFS + n) =  prob_parm.rho_amb * prob_parm.massfrac[n];

  }
  
  //ZND solution 
  else
  {
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real P, rho, U, T, eint, rho_comp, eint_comp;
  
  amrex::Real x_reversed = 1*(x-prob_parm.x_ZND)/100;
  pmf(x_reversed, x_reversed, pmf_vals, prob_parm);
  
  //MKS->CGS conversion
  P = pmf_vals[0]*10; 
  T = pmf_vals[1];
  rho = pmf_vals[2]/1000;
  U = 1.0*pmf_vals[3] * 100; //initialize in shock frame
  
  //assign massfracs
  amrex::Real sum = 0;
  for (int n = 0; n < NUM_SPECIES; n++)
  {
    massfrac[n] = pmf_vals[4+n];
   // sum = sum + massfrac[n];
  }

  //correct the last specie
 // std::cout<<"\nsum = "<<sum; 
 // sum = sum - pmf_vals[4+NUM_SPECIES];
  
  //compute other params
  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2E(T, massfrac, eint);
  eos.PYT2RE(P, massfrac, T, rho_comp, eint_comp);
  
  //assign to state
  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = 1*rho*std::abs(U);
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = rho*eint;
  state(i, j, k, UEDEN) = rho*(eint + 0.5*std::abs(U)*std::abs(U));
  state(i, j, k, UTEMP) = T;

  for (int n = 0; n < NUM_SPECIES; n++){
          if  (massfrac[n]<0)
    {
      std::cout<<"MASSFRAC NEGATIVE!!!!\n";
    }
    state(i, j, k, UFS + n) = rho*massfrac[n];
  }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& /*turb_fluc*/)
{
   if (idir==0 && sgn==1)
  {
    s_ext[URHO]  = prob_parm.rho_amb;
    s_ext[UMX]   = 1.0 * prob_parm.rho_amb*prob_parm.D_CJ;
    s_ext[UMY]   = 0.0;
    s_ext[UMZ]   = 0.0;
    s_ext[UEINT] = prob_parm.rho_amb * prob_parm.eint_amb;
    s_ext[UEDEN] = prob_parm.rho_amb * (prob_parm.eint_amb + 0.5 * prob_parm.D_CJ * prob_parm.D_CJ +
                                                             0.5 * s_ext[UMY] * s_ext[UMY] +
                                                             0.5 * s_ext[UMZ] * s_ext[UMZ]);
    s_ext[UTEMP] = prob_parm.T_amb;
    for (int n = 0; n < NUM_SPECIES; n++){
      s_ext[UFS + n] = prob_parm.rho_amb * prob_parm.massfrac[n]; 
    } 
  }
  else if (idir==1 || idir==2) //move the walls at D_CJ
  {
      s_ext[URHO]  = s_int[URHO];
      s_ext[UMX]   = s_int[URHO]*prob_parm.D_CJ;
      s_ext[UMY]   = -s_int[UMY];
      s_ext[UMZ]   = s_int[UMZ];
      s_ext[UEINT] = s_int[UEINT];

      amrex::Real u, v, w, e_tot_DCJ;
      u = s_ext[UMX]/s_ext[URHO];
      v = s_ext[UMY]/s_ext[URHO];
      w = s_ext[UMZ]/s_ext[URHO];

      e_tot_DCJ = s_int[URHO] * (s_int[UEINT]/s_int[URHO] + 0.5 * prob_parm.D_CJ * prob_parm.D_CJ + 0.5*v*v + 0.5*w*w);

      s_ext[UEDEN] = e_tot_DCJ;

      s_ext[UTEMP] = s_int[UTEMP];
      for (int n = 0; n < NUM_SPECIES; n++){
         s_ext[UFS + n] = s_int[UFS + n];
       }
  }
}

void pc_prob_close();

using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;

#endif
