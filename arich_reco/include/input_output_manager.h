#ifndef input_output_manager_h
#define input_output_manager_h

#include <string>

class input_output_manager
{
public:
  input_output_manager();
  ~input_output_manager();
  bool initialize(int argc, char** argv);
private:
  std::string output_message;
  std::string the_general_mode;
  std::string the_data_input_mode;
  std::string the_data_input_file;
  std::string the_geometry_file;
};

#endif
