#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "input_output_manager.h"

using namespace std;

input_output_manager::input_output_manager()
:output_message("")
{
  cout << "Created the_input_output_manager." << endl;
}

input_output_manager::~input_output_manager()
{
}

bool input_output_manager::initialize(int argc, char** argv)
{
  bool initialization_succeeded = true;
  
  string command_input_file_name = "";
  
  struct stat buffer;
  
  ifstream command_input_file_stream;

  string line;
  string card;
  string string_card_value;
  int int_card_value;
  
  cout << "File initialziation started in the_input_output_manager." << endl;
  
  if (argc == 1)
  {
    output_message += "Too few command line arguments. Usage:\n";
    output_message += argv[0];
    output_message += " command_input_file.json\n";
    initialization_succeeded = false;
  }
  else if (argc > 2)
  {
    output_message += "Too many command line arguments. Usage:\n";
    output_message += argv[0];
    output_message += " command_input_file.json\n";
    
    initialization_succeeded = false;
  }
  else
  {
    command_input_file_name = argv[1];
    
    if (stat(command_input_file_name.c_str(), &buffer) != 0)
    {
      output_message += "The command_input_file does not exist: ";
      output_message += command_input_file_name;
      output_message += "\n";
      
      initialization_succeeded = false;
    }
    else
    {
      output_message += "Using command input file name: ";
      output_message += command_input_file_name;
      output_message += "\n";
      
      command_input_file_stream.open(command_input_file_name.c_str(), ifstream::in);
    
      while (getline(command_input_file_stream, line))
      {
        bool empty_line = false;
        
        if (line.empty())
        {
          empty_line = true;
        }
        else
        {
          bool line_of_spaces = true;
          
          for (int i = 0; i < line.length(); i++)
          {
            if (line[i] != ' ')
            {
              line_of_spaces = false;
            }
          }
          
          if (!line_of_spaces) empty_line = false;
        }
        
        if (!empty_line)
        {
          if (line[0] != '#')
          {
            istringstream temp_string_stream(line);
            temp_string_stream >> card;
            if (card == "general_mode")
            {
              temp_string_stream >> the_general_mode;
              output_message += "Using general mode: ";
              output_message += the_general_mode;
              output_message += "\n";
            }
            else if (card == "data_input_mode")
            {
              temp_string_stream >> the_data_input_mode;
              output_message += "Using data input mode: ";
              output_message += the_data_input_mode;
              output_message += "\n";
            }
            else if (card == "data_input_file")
            {
              temp_string_stream >> the_data_input_file;
              output_message += "Using data input file: ";
              output_message += the_data_input_file;
              output_message += "\n";
            }
            else if (card == "geometry_file")
            {
              temp_string_stream >> the_geometry_file;
              output_message += "Using geomtry file: ";
              output_message += the_geometry_file;
              output_message += "\n";
            }
            else
            {
              output_message += "Unknown input item: ";
              output_message += card;
              output_message += "\n";
              output_message += "Ignored...\n";
            }
          }
        }
      }
    }
    
    command_input_file_stream.close();
  }
  
  cout << output_message;
  
  return initialization_succeeded;
}
