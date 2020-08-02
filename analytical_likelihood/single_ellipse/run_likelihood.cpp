#include <math.h>
#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TApplication.h>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TColor.h>

using namespace std;
using namespace ROOT::Math;

void my_application(int argc, char** argv)
{
  const double PI = 3.14159265359;
  const double to_radians = PI / 180.;
  const double to_degrees = 180. / PI;
  
  const int number_of_points = 12000;
  
  double refraction_angle_1, refraction_angle_2;

  const double refraction_index_1 = 1.0352; // aerogel 1 index of refraction
  const double refraction_index_2 = 1.0452;
  //const double refraction_index_1 = 1.2; // arbitrary large number for verification purposes
  //const double refraction_index_2 = 1.5;
  const double refraction_index_3 = 1.; // air

  const double particle_momentum = 7.0; // in GeV/c
  
  // const double particle_mass = 0.134977; // pion, in GeV/c^2
   const double particle_mass = 0.493677; // kaon, in GeV/c^2
  // const double particle_mass = 0.938272; // proton, in GeV/c^2
  
  const double particle_beta = 1. / (sqrt(1 + particle_mass * particle_mass / (particle_momentum * particle_momentum)));
  
  if (fabs(1./ (refraction_index_1 * particle_beta)) > 1. || fabs(1./ (refraction_index_2 * particle_beta)) > 1.)
  {
    cout << "Particle will not Cerenkov in materials with the given indices of refraction." << endl;
    cout << "cos(theta) are: " << fabs(1./ (refraction_index_1 * particle_beta)) << " " << fabs(1./ (refraction_index_2 * particle_beta)) << endl;
    cout << "Quitting..." << endl;
    exit(2001);
  }
  
  double theta_cone = acos(1./ (refraction_index_1 * particle_beta)); // between axis and radiation
  double cone_direction_angle = 0.2; // 100 milliradians as a typical case
  
  /*
  cout << "particle momentum beta angle" << endl;
  cout << particle_mass << " " << particle_momentum << " " << particle_beta << " " << theta_cone * to_degrees << endl;
  */
  
  const double aerogel_1_thickness = 0.02; // 2 cm
  const double aerogel_2_thickness = 0.02; // 2 cm
  
  const double z_interface_1 = aerogel_1_thickness;
  const double z_interface_2 = aerogel_1_thickness + aerogel_2_thickness;
  const double z_screen = 0.2;
  
  const XYZVector plane_normal(0., 0., 1.);
  
  double rightmost_point_1;
  double leftmost_point_1;
  double upmost_point_1;

  XYZVector cone_vectors[number_of_points];
  XYZVector hit_points_1[number_of_points];
  XYZVector refracted_cone_vectors_1[number_of_points];
  XYZVector hit_points_2[number_of_points];
  XYZVector refracted_cone_vectors_2[number_of_points];
  XYZVector hit_points_3[number_of_points];
  
  RotationY ry(cone_direction_angle);

  const double coef_R_squared = tan(theta_cone) * tan(theta_cone);
  const double coef_A = sin(cone_direction_angle);
  const double coef_C = - cos(cone_direction_angle);
  const double coef_D = aerogel_1_thickness;
  
  double max_x_1, min_x_1, max_x_2, min_x_2, max_x_3, min_x_3;
  double max_y_1, min_y_1, max_y_2, min_y_2, max_y_3, min_y_3;
  int id_max_1, id_min_1, id_max_2, id_min_2, id_max_3, id_min_3; // these are for y only
  
  TMultiGraph *mg = new TMultiGraph();
  TCanvas *c1 = new TCanvas("c1", "c1", 5, 26, 1200, 815);
  TGraph *gr11 = new TGraph(number_of_points);
  TGraph *gr12 = new TGraph(number_of_points);
  TGraph *gr13 = new TGraph(number_of_points);
  
  TGraph *gr21 = new TGraph(3);

  double dummy_double;
  
  ofstream out_file("output.txt");
  
  for (int i = 0; i < number_of_points; i++)
  {
    // cone vectors distributed uniformly
    dummy_double = ((double) i) * 2. * PI / ((double) number_of_points);
    
    // set the cone direction lines as unit vectors
    cone_vectors[i].SetCoordinates(sin(theta_cone) * cos(dummy_double), sin(theta_cone) * sin(dummy_double), cos(theta_cone));
    
    cone_vectors[i] = ry * cone_vectors[i];
    
    // project on the first plane
    hit_points_1[i] = cone_vectors[i] * z_interface_1 / fabs(cone_vectors[i].Dot(plane_normal));
    
    
    // cone_vectors are still unit
    // cout << cone_vectors[i].Mag2() << endl;
    refraction_angle_1 = acos(fabs(cone_vectors[i].Dot(plane_normal)));
    refraction_angle_2 = asin((refraction_index_1 / refraction_index_2) * sin(refraction_angle_1));
    refracted_cone_vectors_1[i].SetCoordinates(cone_vectors[i].x(), cone_vectors[i].y(), 0);
    refracted_cone_vectors_1[i] = refracted_cone_vectors_1[i].Unit();
    refracted_cone_vectors_1[i] = refracted_cone_vectors_1[i] * sin(refraction_angle_2);
    refracted_cone_vectors_1[i].SetZ(cos(refraction_angle_2));
    
    hit_points_2[i] = hit_points_1[i] + refracted_cone_vectors_1[i] * (z_interface_2 - z_interface_1) / fabs(refracted_cone_vectors_1[i].Dot(plane_normal));
    
    refraction_angle_1 = acos(fabs(refracted_cone_vectors_1[i].Dot(plane_normal)));
    refraction_angle_2 = asin((refraction_index_2 / refraction_index_3) * sin(refraction_angle_1));
    refracted_cone_vectors_2[i].SetCoordinates(refracted_cone_vectors_1[i].x(), refracted_cone_vectors_1[i].y(), 0);
    refracted_cone_vectors_2[i] = refracted_cone_vectors_2[i].Unit();
    refracted_cone_vectors_2[i] = refracted_cone_vectors_2[i] * sin(refraction_angle_2);
    refracted_cone_vectors_2[i].SetZ(cos(refraction_angle_2));
    
    hit_points_3[i] = hit_points_2[i] + refracted_cone_vectors_2[i] * (z_screen - z_interface_2) / fabs(refracted_cone_vectors_2[i].Dot(plane_normal));
    
    // so hit_points_1 are now the positions of the hits in the lab system (particle starts at (0, 0, 0))
    gr11->SetPoint(i, hit_points_1[i].x(), hit_points_1[i].y());
    
    gr12->SetPoint(i, hit_points_2[i].x(), hit_points_2[i].y());
    
    gr13->SetPoint(i, hit_points_3[i].x(), hit_points_3[i].y());
  }
  
  leftmost_point_1  = aerogel_1_thickness * tan(cone_direction_angle - theta_cone);
  rightmost_point_1 = aerogel_1_thickness * tan(cone_direction_angle + theta_cone);
  
  gr21->SetPoint(0, leftmost_point_1, 0.);
  gr21->SetPoint(1, (leftmost_point_1 + rightmost_point_1) / 2., 0.);
  gr21->SetPoint(2, rightmost_point_1, 0.);
  
  Int_t ci_11 = TColor::GetFreeColorIndex();
  TColor* the_color_11 = new TColor(ci_11, 0., 0., 139./255.);
  Int_t ci_12 = TColor::GetFreeColorIndex();
  TColor* the_color_12 = new TColor(ci_12, 0., 0., 1.);
  Int_t ci_13 = TColor::GetFreeColorIndex();
  TColor* the_color_13 = new TColor(ci_13, 0., 191./255., 1.);
  Int_t ci_21 = TColor::GetFreeColorIndex();
  TColor* the_color_21 = new TColor(ci_21, 1., 0., 0.);
  
  dummy_double = (coef_R_squared * coef_A * coef_D) / (coef_C * coef_C - coef_R_squared * coef_A * coef_A);
  upmost_point_1 = sqrt((coef_R_squared * coef_A * coef_A - coef_C * coef_C) * dummy_double * dummy_double +
                        2. * coef_R_squared * coef_A * coef_D * dummy_double +
                        coef_R_squared * coef_D * coef_D) / coef_C;
  
  /*
  cout << upmost_point_1 << endl;
  cout << hit_points_1[number_of_points / 4 - 1].y() << endl;
  
  cout << "Semimajor axis: " << (rightmost_point_1 - leftmost_point_1) / 2. << endl;
  */
  
  // the initialization is assuming the highest y is positive and the lowest y is negative
  // this is due to the choice of the local coordinate system
  max_y_1 = 0.;
  min_y_1 = 0.;
  max_y_2 = 0.;
  min_y_2 = 0.;
  max_y_3 = 0.;
  min_y_3 = 0.;
  
  max_x_1 = -1000.;
  min_x_1 = 1000.;
  max_x_2 = -1000.;
  min_x_2 = 1000.;
  max_x_3 = -1000.;
  min_x_3 = 1000.;
  
  for (int i = 0; i < number_of_points; i++)
  {
    //cout << i << " " << hit_points_1[i].y() << " " << hit_points_2[i].y() << endl;
    
    if (hit_points_1[i].y() > max_y_1)
    {
      max_y_1 = hit_points_1[i].y();
      id_max_1 = i;
    }
    
    if (hit_points_1[i].y() < min_y_1)
    {
      min_y_1 = hit_points_1[i].y();
      id_min_1 = i;
    }
    
    if (hit_points_2[i].y() > max_y_2)
    {
      max_y_2 = hit_points_2[i].y();
      id_max_2 = i;
    }
    
    if (hit_points_2[i].y() < min_y_2)
    {
      min_y_2 = hit_points_2[i].y();
      id_min_2 = i;
    }
    
    if (hit_points_3[i].y() > max_y_3)
    {
      max_y_3 = hit_points_3[i].y();
      id_max_3 = i;
    }
    
    if (hit_points_3[i].y() < min_y_3)
    {
      min_y_3 = hit_points_3[i].y();
      id_min_3 = i;
    }
    
    if (hit_points_1[i].x() > max_x_1)
    {
      max_x_1 = hit_points_1[i].x();
    }
    
    if (hit_points_1[i].x() < min_x_1)
    {
      min_x_1 = hit_points_1[i].x();
    }
    
    if (hit_points_2[i].x() > max_x_2)
    {
      max_x_2 = hit_points_2[i].x();
    }
    
    if (hit_points_2[i].x() < min_x_2)
    {
      min_x_2 = hit_points_2[i].x();
    }
    
    if (hit_points_3[i].x() > max_x_3)
    {
      max_x_3 = hit_points_3[i].x();
    }
    
    if (hit_points_3[i].x() < min_x_3)
    {
      min_x_3 = hit_points_3[i].x();
    }
  }

  for (int i = 0; i < number_of_points; i++)
  {
    out_file << hit_points_3[i].x() * 100. << " " << hit_points_3[i].y() * 100. << endl;
  }
  
  out_file.close();
  
  cout << "Maximum y index: " << id_max_1 << " " << id_max_2 << " " << id_max_3 << endl;
  cout << "Minimum y index: " << id_min_1 << " " << id_min_2 << " " << id_min_3 << endl;

  cout << "Eccentricities: " << sqrt(1 - ((max_y_1 - min_y_1) * (max_y_1 - min_y_1)) / ((max_x_1 - min_x_1) * (max_x_1 - min_x_1))) << " "
                             << sqrt(1 - ((max_y_2 - min_y_2) * (max_y_2 - min_y_2)) / ((max_x_2 - min_x_2) * (max_x_2 - min_x_2))) << " "
                             << sqrt(1 - ((max_y_3 - min_y_3) * (max_y_3 - min_y_3)) / ((max_x_3 - min_x_3) * (max_x_3 - min_x_3))) << endl;
  
  gr21->SetPoint(3, (leftmost_point_1 + rightmost_point_1) / 2., upmost_point_1);
  gr21->SetPoint(4, (leftmost_point_1 + rightmost_point_1) / 2., -upmost_point_1);
  
  gr21->SetPoint(5, z_screen * tan(cone_direction_angle), 0.);

  gr11->SetMarkerColor(ci_11);
  gr11->SetMarkerStyle(20);
  gr11->SetMarkerSize(0.2);
  
  gr12->SetMarkerColor(ci_12);
  gr12->SetMarkerStyle(20);
  gr12->SetMarkerSize(0.2);
  
  gr13->SetMarkerColor(ci_13);
  gr13->SetMarkerStyle(20);
  gr13->SetMarkerSize(0.2);
  
  gr21->SetMarkerColor(ci_21);
  gr21->SetMarkerStyle(20);
  gr21->SetMarkerSize(0.6);
  
  mg->Add(gr11);
  mg->Add(gr12);
  mg->Add(gr13);
  mg->Add(gr21);

  c1->cd();
  mg->Draw("AP");
  
  return;
}

int main(int argc, char** argv)
{
  TApplication app("EMPHATIC ARICH Analytical Reconstruction", &argc, argv);
  my_application(app.Argc(), app.Argv());
  app.Run();
  
  // app.Terminate(); // still doesn't quit
  
  return 0;
}



