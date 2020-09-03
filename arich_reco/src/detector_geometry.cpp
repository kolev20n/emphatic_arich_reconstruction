#include "detector_geometry.h"

detector_geometry::detector_geometry(int array_rows, int array_columns)
: number_of_aerogel_layers(2),
  number_of_detector_array_rows(array_rows),
  number_of_detector_array_columns(array_columns)
{
}

detector_geometry::~detector_geometry()
{
}


const int detector_geometry::get_number_of_aerogel_layers() const
{
  return number_of_aerogel_layers;
}

const int detector_geometry::get_number_of_detector_array_rows() const
{
  return number_of_detector_array_rows;
}

const int detector_geometry::get_number_of_detector_array_columns() const
{
  return number_of_detector_array_columns;
}

