/* Copyright 2015 Dominic Deuber
*
* This file is part of VCFvifia.
*
* VCFvifia is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* VCFvifia is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with VCFvifia. If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file utils.hpp
 * collection of useful functions
 */

#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

namespace utils{
  /**
   * function for printing the elements of a vector to a file
   * @param aVector the vector to be printed
   * @param aFilename the name of the output file
   * @param aHeaderline optional header line to be printed
   * @return true in case of success, false otherwise
   */
  template<typename T>
  inline bool printVectorToFile(std::vector<T> &aVector, const std::string &aFilename, const std::string &aHeaderline="") {
    if(!aVector.empty()) {
      std::ofstream outfile(aFilename);
      if(!outfile.is_open()) {
        std::cerr << "Warning: could not create file: '" << aFilename << "'" << std::endl;
        return false;
      }
      if(aHeaderline != "") {
        outfile << aHeaderline << std::endl;
      }
      for(auto itr = aVector.begin(); itr != aVector.end(); ++itr) {
        outfile << *itr << std::endl;
      }
      outfile.close();
    }
    return true;
  }

  /**
   * function for getting the string representation of a vector
   * @param aVector the vector
   * @return string representation in case of success, empty string otherwise
   */
  template<typename T>
  inline std::string toString(std::vector<T> &aVector) {
    std::stringstream sstream;
    for(auto itr = aVector.begin(); itr != aVector.end(); ++itr) {
      sstream << (*itr) << " ";
    }
    return sstream.str();
  }

  /**
   * function for checking wheter first char of a given string is digit
   * @param aString the string to be checked
   * @return true in case of first char being digit, false otherwise
   */
  bool firstIsDigit(const std::string & aString);
}
#endif
