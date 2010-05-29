#ifndef GRAPHIC_PARAMETERS_H
#define GRAPHIC_PARAMETERS_H


class graphic_parameters
{ public:
int show_pointcharges;
double scale_pointcharges;
double show_abc_unitcell,show_primitive_crystal_unitcell,show_magnetic_unitcell;
double show_atoms,scale_view_1,scale_view_2,scale_view_3;
double show_chargedensity,spins_scale_moment;
double spins_wave_amplitude;
int showprim;
double spins_show_ellipses;
double spins_scale_static_moment;
double spins_show_static_moment_direction;
double spins_show_oscillation;
double threshhold;
char * title;
~graphic_parameters(){delete title;}
graphic_parameters()
{title=new char[100];
 sprintf(title,"chargedensity");
show_pointcharges=0;scale_pointcharges=1;
showprim=0;
show_abc_unitcell=1.0;
show_primitive_crystal_unitcell=1.0;
show_magnetic_unitcell=1.0;
show_atoms=1.0;
scale_view_1=1.0;
scale_view_2=1.0;
scale_view_3=1.0;
show_chargedensity=1.0;
spins_scale_moment=1.0;
spins_show_ellipses=1;
spins_scale_static_moment=1;
spins_show_static_moment_direction=1;
spins_wave_amplitude=0;
spins_show_oscillation=0;
threshhold=0.05;
}

};
#endif