#include <iostream>

#include "input_output_manager.h"
#include "detector_geometry.h"

using namespace std;

int main(int argc, char** argv)
{
  input_output_manager* the_input_output_manager = new input_output_manager();
  
  if (!the_input_output_manager->initialize(argc, argv))
  {
    cout << "The initialize() of the_input_output_manager failed. Quitting..." << endl;
    
    return(101);
  }
  
  int a = 3;
  int b = 3;
  
  detector_geometry* the_detector_geometry = new detector_geometry(a, b);
  cout << the_detector_geometry->get_number_of_aerogel_layers() << endl;
  cout << the_detector_geometry->get_number_of_detector_array_rows() << endl;
  cout << the_detector_geometry->get_number_of_detector_array_columns() << endl;

  return 0;
}
