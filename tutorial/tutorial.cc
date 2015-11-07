/*
 * This is the pseudo-code for the tutorial step to implement in deal.II.
 */


void setup()
{
  // I think what would be cool if we could do a projection on a surface: we can
  // use hyper_cube_with_cylindrical_hole and a cylinder. This is what you want
  // to do when you study a fuel assembly.
  // We should make it so that the meshes on the surfaces do not share the same
  // cells.

  GridGenerator::hyper_cube_with_cylindrical_hole(fluid_tria,radius,side,
      length,repetition);
  global_refine(fluid_tria,4);
  GridGenerator::cylinder(solid_tria,radius,length/2.);
  // Shift the cylinder
  Tensor<1,3> shift_vector;
  shift_vector[2] = length/2.;
  GridTools::shift(shift_vector,solid_tria);
  global_refine(solid_tria,4);

  // Set a function on solid_tria that we will project on fluid_tria. Something
  // like f(z) = a*z;
  Vector<double> solid_temperature = compute_solid_temperature();
}


void project()
{
  // Inputs: solid_dof_handler, solid_temperature, and liquid_dof_handler.
  // Output: liquid_temperature
  project(solid_dof_handler,solid_temperature,liquid_dof_handler,
      liquid_temperature);
}

int main()
{
  setup();
  
  project();

  output();

  return 0;
}
