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
 * @file Plotter.hpp
 * functions to plot data
 */

#ifndef PLOTTER_HPP
#define PLOTTER_HPP

#include "gnuplot-iostream.h"
#include <boost/filesystem.hpp>

/**
 * plotting class
 */
class CPlotter {
public:
  /**
   * struct holding basic format options of the plots
   */
  struct FormatSpec {
    std::string title; /**< title of the plot */
    std::string xlabel; /**< label of the x-axis */
    std::string ylabel; /**< label of the y-axis */
    std::string y2label; /**< label of the right y-axis */
    std::string suffix; /**< suffix for the filename */
    FormatSpec(std::string aTitle, std::string aXLabel, std::string aYLabel, std::string aYLabel2, std::string aSuffix) : title(aTitle), xlabel(aXLabel), ylabel(aYLabel), y2label(aYLabel2), suffix(aSuffix) { }
  };

  /**
   * constructor of the plotting class
   * @param aFileformat specify the format the plots will be saved in
   */
  CPlotter(const std::string &aFileformat);

  /**
   * function to plot a bar graph
   * @param aData data to be plotted (if empty, no plot will be generated)
   * @param aFilename name to create a filename for the resulting plot
   * @param aFormat struct to specify the format of the plots
   */
  template<typename X, typename Y> void plotBarGraph(std::vector<std::pair<X, Y>> aData, const std::string &aFilename, FormatSpec aFormat) {
    if(aData.empty()) {
      std::cout << "Warning: no data to create plot: " << aFormat.title << std::endl;
    } else {
      Gnuplot gp;
      preSettings(gp, aFilename, aFormat, "bmargin");
      gp << "set style fill solid\n";
      gp << "set xtic 1\n";
      gp << "set boxwidth 0.75\n";
      gp << "set nokey\n";
      gp << "plot '-' using 2:xtic(1) with boxes lt rgb \"blue\"\n";
      gp.send1d(aData);
    }
  }

  /**
   * function to plot a distribution graph
   * @param aData data to be plotted (if empty, no plot will be generated)
   * @param aFilename name to create a filename for the resulting plot
   * @param aFormat struct to specify the format of the plots
   */
  template<typename X, typename Y>void plotDistributionGraph(std::vector<std::pair<X, Y>> aData, const std::string &aFilename, FormatSpec aFormat) {
    if(aData.empty()) {
      std::cout << "Warning: no data to create plot: " << aFormat.title << std::endl;
    } else {
      Gnuplot gp;
      preSettings(gp, aFilename, aFormat, "bmargin");
      gp << "set style fill solid\n";
      gp << "set xtic 4\n";
      gp << "set xtics rotate by -45\n";
      gp << "set boxwidth 0.25\n";
      gp << "set nokey\n";
      gp << "plot '-' using 1:2 smooth frequency with boxes lt rgb \"blue\"\n";
      gp.send1d(aData);
    }
  }

  /**
   * function to plot a stacked bar graph
   * @param aDataPath path to the data to be plotted (if file is empty, no plots were generated)
   * @param aFilename name to create a filename for the resulting plot
   * @param aFormat struct to specify the format of the plots
   */
  void plotStackedBarGraph(const std::string& aDataPath, const std::string &aFilename, FormatSpec aFormat);

  /**
   * function to plot a bar graph with two bars
   * @param aDataPath path to the data to be plotted (if file is empty, no plots were generated)
   * @param aFilename name to create a filename for the resulting plot
   * @param aFormat struct to specify the format of the plots
   */
  void plotTwoBarsGraph(const std::string& aDataPath, const std::string &aFilename, FormatSpec aFormat);

  /**
   * function to plot a weighted distribution bar graph
   * @param aDataPath path to the data to be plotted (if file is empty, no plots were generated)
   * @param aFilename name to create a filename for the resulting plot
   * @param aFormat struct to specify the format of the plots
   */
  void plotWeightedDistributionGraph(const std::string& aDataPath, const std::string &aFilename, FormatSpec aFormat);

private:
  /*
   * fileformat the plots will be saved in
   */
  std::string mFileformat;

  /*
   * function to set basic information to gnuplot
   * @param aGnuplot Gnuplot object where to set the information
   * @param aFilename name to create a filename for the resulting plot
   * @param aFormat struct to specify the format of the plots
   * @param aKeyposition string to define the position of the legend
   */
  void preSettings(Gnuplot &aGnuplot, const std::string &aFilename, FormatSpec aFormat, std::string aKeyposition);
};
#endif //PLOTTER_HPP
