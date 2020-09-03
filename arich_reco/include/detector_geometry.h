#ifndef detector_geometry_h
#define detector_geometry_h

#include <vector>

class detector_geometry
{
public:
  detector_geometry(int array_rows, int array_columns);
  ~detector_geometry();
  const int get_number_of_aerogel_layers() const;
  const int get_number_of_detector_array_rows() const;
  const int get_number_of_detector_array_columns() const;
private:
  const int number_of_aerogel_layers;
  const int number_of_detector_array_rows;
  const int number_of_detector_array_columns;
  std::vector<std::vector<double>> aerogel_size;
  // material ??
  //individual_detector_cell the_detector_array
  std::vector<std::vector<int>> the_detector_array;
};

#endif
