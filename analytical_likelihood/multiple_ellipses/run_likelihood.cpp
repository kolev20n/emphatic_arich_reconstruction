#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

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
  
  const int number_of_points = 1200;
  const int number_of_layers = 2; // aerogel layers
  
  const double refraction_index_aerogel[number_of_layers] = {1.0352, 1.0452}; // aerogel indices of refraction
  const double refraction_index_air = 1.;
 
  const double aerogel_thickness[number_of_layers] = {0.02, 0.02}; // 1 cm each layer
  
  const double z_interface[number_of_layers] = {aerogel_thickness[0], aerogel_thickness[0] + aerogel_thickness[1]};
  
  const XYZVector plane_normal(0., 0., 1.);
  
  XYZVector cone_vectors_aerogel_1[number_of_points];
  XYZVector cone_vectors_aerogel_1_refracted_once[number_of_points];
  XYZVector cone_vectors_aerogel_1_refracted_twice[number_of_points];

  XYZVector cone_vectors_aerogel_2[number_of_points];
  XYZVector cone_vectors_aerogel_2_refracted_once[number_of_points];
  
  XYZVector hit_points_aerogel_1_front_on_interface_1[number_of_points];
  XYZVector hit_points_aerogel_1_front_on_interface_2[number_of_points];
  XYZVector hit_points_aerogel_1_front_on_screen[number_of_points];
  XYZVector hit_points_aerogel_1_back_on_interface_2[number_of_points];
  XYZVector hit_points_aerogel_1_back_on_screen[number_of_points];
  
  XYZVector hit_points_aerogel_2_front_on_interface_2[number_of_points];
  XYZVector hit_points_aerogel_2_front_on_screen[number_of_points];
  XYZVector hit_points_aerogel_2_back_on_screen[number_of_points];
  
  double particle_momentum = 5.0; // in GeV/c; default value
  double particle_mass = 0.493677; // kaon, in GeV/c^2; default value
  // const double particle_mass = 0.134977; // pion, in GeV/c^2
  // const double particle_mass = 0.938272; // proton, in GeV/c^2
  double cone_direction_angle = 0.1; // 100 milliradians as a typical case
  double z_screen = 0.2; // in m
  
  string a_card;
  double a_value;
  
  ifstream in_file(argv[1]);
  
  while (in_file >> a_card >> a_value)
  {
    if (a_card == "particle_momentum")
    {
      particle_momentum = a_value;
    }
    else if (a_card == "particle_mass")
    {
      particle_mass = a_value;
    }
    else if (a_card == "cone_direction_angle")
    {
      cone_direction_angle = a_value;
    }
    else if (a_card == "z_screen")
    {
      z_screen = a_value;
    }
  }
  
  RotationY ry(cone_direction_angle);

  double particle_beta = 1. / (sqrt(1. + particle_mass * particle_mass / (particle_momentum * particle_momentum)));
  
  double refraction_angle[number_of_layers];
  
  double theta_cone[number_of_layers]; // between axis and radiation
  theta_cone[0] = acos(1./ (refraction_index_aerogel[0] * particle_beta));
  theta_cone[1] = acos(1./ (refraction_index_aerogel[1] * particle_beta));

  XYZVector particle_hit_point_on_interface_1, particle_hit_point_on_interface_2;
  particle_hit_point_on_interface_1.SetCoordinates(aerogel_thickness[0] * sin(cone_direction_angle), 0., 0.);
  particle_hit_point_on_interface_2.SetCoordinates((aerogel_thickness[0] + aerogel_thickness[1]) * sin(cone_direction_angle), 0., 0.);

  TMultiGraph *mg = new TMultiGraph();
  
  TCanvas *c1 = new TCanvas("c1", "c1", 5, 26, 1200, 815);
  TGraph *gr11 = new TGraph(number_of_points);
  TGraph *gr12 = new TGraph(number_of_points);
  TGraph *gr13 = new TGraph(number_of_points);
  
  TGraph *gr22 = new TGraph(number_of_points);
  TGraph *gr23 = new TGraph(number_of_points);
  
  TGraph *gr32 = new TGraph(number_of_points);
  TGraph *gr33 = new TGraph(number_of_points);
  
  TGraph *gr43 = new TGraph(number_of_points);
  
  double dummy_double;
  
  for (int i = 0; i < number_of_points; i++)
  {
    // cone vectors distributed uniformly
    dummy_double = ((double) i) * 2. * PI / ((double) number_of_points);
    
    // set the cone direction lines as unit vectors
    cone_vectors_aerogel_1[i].SetCoordinates(sin(theta_cone[0]) * cos(dummy_double), sin(theta_cone[0]) * sin(dummy_double), cos(theta_cone[0]));
    cone_vectors_aerogel_2[i].SetCoordinates(sin(theta_cone[1]) * cos(dummy_double), sin(theta_cone[1]) * sin(dummy_double), cos(theta_cone[1]));
    
    cone_vectors_aerogel_1[i] = ry * cone_vectors_aerogel_1[i];
    cone_vectors_aerogel_2[i] = ry * cone_vectors_aerogel_2[i];

    // project on the first plane
    hit_points_aerogel_1_front_on_interface_1[i] = cone_vectors_aerogel_1[i] * z_interface[0] / fabs(cone_vectors_aerogel_1[i].Dot(plane_normal));
    
    // cone_vectors are still unit
    // cout << cone_vectors[i].Mag2() << endl;
    refraction_angle[0] = acos(fabs(cone_vectors_aerogel_1[i].Dot(plane_normal)));
    refraction_angle[1] = asin((refraction_index_aerogel[0] / refraction_index_aerogel[1]) * sin(refraction_angle[0]));
    cone_vectors_aerogel_1_refracted_once[i].SetCoordinates(cone_vectors_aerogel_1[i].x(), cone_vectors_aerogel_1[i].y(), 0.);
    cone_vectors_aerogel_1_refracted_once[i] = cone_vectors_aerogel_1_refracted_once[i].Unit();
    cone_vectors_aerogel_1_refracted_once[i] = cone_vectors_aerogel_1_refracted_once[i] * sin(refraction_angle[1]);
    cone_vectors_aerogel_1_refracted_once[i].SetZ(cos(refraction_angle[1]));
    
    hit_points_aerogel_1_front_on_interface_2[i] = hit_points_aerogel_1_front_on_interface_1[i] + cone_vectors_aerogel_1_refracted_once[i] * (z_interface[1] - z_interface[0]) / fabs(cone_vectors_aerogel_1_refracted_once[i].Dot(plane_normal));
    
    hit_points_aerogel_1_back_on_interface_2[i] = particle_hit_point_on_interface_1 + cone_vectors_aerogel_1_refracted_once[i] * (z_interface[1] - z_interface[0]) / fabs(cone_vectors_aerogel_1_refracted_once[i].Dot(plane_normal));
    
    refraction_angle[0] = acos(fabs(cone_vectors_aerogel_1_refracted_once[i].Dot(plane_normal)));
    refraction_angle[1] = asin((refraction_index_aerogel[1] / refraction_index_air) * sin(refraction_angle[0]));
    cone_vectors_aerogel_1_refracted_twice[i].SetCoordinates(cone_vectors_aerogel_1_refracted_once[i].x(), cone_vectors_aerogel_1_refracted_once[i].y(), 0.);
    cone_vectors_aerogel_1_refracted_twice[i] = cone_vectors_aerogel_1_refracted_twice[i].Unit();
    cone_vectors_aerogel_1_refracted_twice[i] = cone_vectors_aerogel_1_refracted_twice[i] * sin(refraction_angle[1]);
    cone_vectors_aerogel_1_refracted_twice[i].SetZ(cos(refraction_angle[1]));
    
    refraction_angle[0] = acos(fabs(cone_vectors_aerogel_2[i].Dot(plane_normal)));
    refraction_angle[1] = asin((refraction_index_aerogel[1] / refraction_index_air) * sin(refraction_angle[0]));
    cone_vectors_aerogel_2_refracted_once[i].SetCoordinates(cone_vectors_aerogel_2[i].x(), cone_vectors_aerogel_2[i].y(), 0.);
    cone_vectors_aerogel_2_refracted_once[i] = cone_vectors_aerogel_2_refracted_once[i].Unit();
    cone_vectors_aerogel_2_refracted_once[i] = cone_vectors_aerogel_2_refracted_once[i] * sin(refraction_angle[1]);
    cone_vectors_aerogel_2_refracted_once[i].SetZ(cos(refraction_angle[1]));
    
    hit_points_aerogel_1_front_on_screen[i] = hit_points_aerogel_1_front_on_interface_2[i] + cone_vectors_aerogel_1_refracted_twice[i] * (z_screen - z_interface[1]) / fabs(cone_vectors_aerogel_1_refracted_twice[i].Dot(plane_normal));
    
    hit_points_aerogel_1_back_on_screen[i] = hit_points_aerogel_1_back_on_interface_2[i] + cone_vectors_aerogel_1_refracted_twice[i] * (z_screen - z_interface[1]) / fabs(cone_vectors_aerogel_1_refracted_twice[i].Dot(plane_normal));
    
    hit_points_aerogel_2_front_on_interface_2[i] = particle_hit_point_on_interface_1 + cone_vectors_aerogel_2[i] * (z_interface[1] - z_interface[0]) / fabs(cone_vectors_aerogel_2[i].Dot(plane_normal));
    
    hit_points_aerogel_2_front_on_screen[i] = hit_points_aerogel_2_front_on_interface_2[i] + cone_vectors_aerogel_2_refracted_once[i] * (z_screen - z_interface[1]) / fabs(cone_vectors_aerogel_2_refracted_once[i].Dot(plane_normal));
    
    hit_points_aerogel_2_back_on_screen[i] = particle_hit_point_on_interface_2 + cone_vectors_aerogel_2_refracted_once[i] * (z_screen - z_interface[1]) / fabs(cone_vectors_aerogel_2_refracted_once[i].Dot(plane_normal));
    
    // so hit_points_1 are now the positions of the hits in the lab system (particle starts at (0, 0, 0))
    gr11->SetPoint(i, hit_points_aerogel_1_front_on_interface_1[i].x(), hit_points_aerogel_1_front_on_interface_1[i].y());
    
    gr12->SetPoint(i, hit_points_aerogel_1_front_on_interface_2[i].x(), hit_points_aerogel_1_front_on_interface_2[i].y());
    
    gr13->SetPoint(i, hit_points_aerogel_1_front_on_screen[i].x(), hit_points_aerogel_1_front_on_screen[i].y());
    
    gr22->SetPoint(i, hit_points_aerogel_1_back_on_interface_2[i].x(), hit_points_aerogel_1_back_on_interface_2[i].y());
    
    gr23->SetPoint(i, hit_points_aerogel_1_back_on_screen[i].x(), hit_points_aerogel_1_back_on_screen[i].y());
    
    gr32->SetPoint(i, hit_points_aerogel_2_front_on_interface_2[i].x(), hit_points_aerogel_2_front_on_interface_2[i].y());
    
    gr33->SetPoint(i, hit_points_aerogel_2_front_on_screen[i].x(), hit_points_aerogel_2_front_on_screen[i].y());
    
    gr43->SetPoint(i, hit_points_aerogel_2_back_on_screen[i].x(), hit_points_aerogel_2_back_on_screen[i].y());
  }
  
  ofstream out_file("kaon_7gevc_0_0_1.txt");
  
  for (int i = 0; i< number_of_points; i++)
  {
    out_file << hit_points_aerogel_1_front_on_screen[i].x() * 100. << " " << hit_points_aerogel_1_front_on_screen[i].y() * 100. << endl;
    out_file << hit_points_aerogel_1_back_on_screen[i].x() * 100. << " " << hit_points_aerogel_1_back_on_screen[i].y() * 100. << endl;
    out_file << hit_points_aerogel_2_front_on_screen[i].x() * 100. << " " << hit_points_aerogel_2_front_on_screen[i].y() * 100. << endl;
    out_file << hit_points_aerogel_2_back_on_screen[i].x() * 100. << " " << hit_points_aerogel_2_back_on_screen[i].y() * 100. << endl;
  }
  
  out_file.close();
  
  Int_t ci_11 = TColor::GetFreeColorIndex();
  TColor* the_color_11 = new TColor(ci_11, 0., 0., 139./255.);
  Int_t ci_12 = TColor::GetFreeColorIndex();
  TColor* the_color_12 = new TColor(ci_12, 0., 0., 1.);
  Int_t ci_13 = TColor::GetFreeColorIndex();
  TColor* the_color_13 = new TColor(ci_13, 0., 191./255., 1.);
  
  Int_t ci_22 = TColor::GetFreeColorIndex();
  TColor* the_color_22 = new TColor(ci_22, 102./255., 0., 204./255.);
  Int_t ci_23 = TColor::GetFreeColorIndex();
  TColor* the_color_23 = new TColor(ci_23, 178./255., 102./255., 1.);

  Int_t ci_32 = TColor::GetFreeColorIndex();
  TColor* the_color_32 = new TColor(ci_32, 0., 153./255., 0.);
  Int_t ci_33 = TColor::GetFreeColorIndex();
  TColor* the_color_33 = new TColor(ci_33, 51./255., 1., 51./255.);
  
  Int_t ci_43 = TColor::GetFreeColorIndex();
  TColor* the_color_43 = new TColor(ci_43, 1., 1., 51./255.);
  
  gr11->SetMarkerColor(ci_11);
  gr11->SetMarkerStyle(20);
  gr11->SetMarkerSize(0.2);
  
  gr12->SetMarkerColor(ci_12);
  gr12->SetMarkerStyle(20);
  gr12->SetMarkerSize(0.2);
  
  gr13->SetMarkerColor(ci_13);
  gr13->SetMarkerStyle(20);
  gr13->SetMarkerSize(0.2);
 
  gr22->SetMarkerColor(ci_22);
  gr22->SetMarkerStyle(20);
  gr22->SetMarkerSize(0.2);
  
  gr23->SetMarkerColor(ci_23);
  gr23->SetMarkerStyle(20);
  gr23->SetMarkerSize(0.2);
  
  gr32->SetMarkerColor(ci_32);
  gr32->SetMarkerStyle(20);
  gr32->SetMarkerSize(0.2);
  
  gr33->SetMarkerColor(ci_33);
  gr33->SetMarkerStyle(20);
  gr33->SetMarkerSize(0.2);
  
  gr43->SetMarkerColor(ci_43);
  gr43->SetMarkerStyle(20);
  gr43->SetMarkerSize(0.2);
  
  mg->Add(gr11);
  mg->Add(gr12);
  mg->Add(gr13);

  mg->Add(gr22);
  mg->Add(gr23);
  
  mg->Add(gr32);
  mg->Add(gr33);
  
  mg->Add(gr43);
  
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



