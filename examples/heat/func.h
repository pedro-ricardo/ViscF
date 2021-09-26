#ifndef FUNC_H_
#define FUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "kernel.h"

void main_main ();

void advance (amrex::MultiFab& phi_old,
              amrex::MultiFab& phi_new,
	      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& flux,
	      amrex::Real dt,
              amrex::Geometry const& geom);

void init_phi (amrex::MultiFab& phi_new, amrex::Geometry const& geom);
#endif