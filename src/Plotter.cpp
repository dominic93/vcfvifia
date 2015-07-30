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

#include "Plotter.hpp"

using namespace std;

CPlotter::CPlotter(const string &aFileformat) {
  mFileformat = aFileformat;
}

void CPlotter::preSettings(Gnuplot &aGnuplot, const string &aFilename, FormatSpec aFormat, string aKeyposition) {
  string filename = aFilename.substr(aFilename.find_last_of("\\/")+1); //remove the path information (in case of aFilename including it)
	aGnuplot << "set key " << aKeyposition << " autotitle columnheader\n";
  aGnuplot << "set terminal " << mFileformat << "\n";
	std::cout << "Created plot: " << aFormat.title << ", see "  << filename << "_" << aFormat.suffix << "." << mFileformat << std::endl;
	aGnuplot << "set output '" << filename << "_" << aFormat.suffix << "." << mFileformat << "'\n";
  aGnuplot << "set title \"" << aFormat.title << "\"\n";
  aGnuplot << "set ylabel \"" << aFormat.ylabel << "\"\n";
  aGnuplot << "set y2label \"" << aFormat.y2label << "\"\n";
  aGnuplot << "set xlabel \"" << aFormat.xlabel << "\"\n";
}

void CPlotter::plotStackedBarGraph(const string &aDataPath, const string &aFilename, FormatSpec aFormat) {
  if ( !boost::filesystem::exists(aDataPath)) {
    cout << "Warning: no data to create plot: " << aFormat.title << endl;
  } else {
    Gnuplot gp;
    preSettings(gp, aFilename, aFormat, "right");
    gp << "set yrange [0:100]\n";
    gp << "set ylabel \"% of total\"\n";
    gp << "set style data histogram\n";
    gp << "set style histogram rowstacked\n";
    gp << "set style fill solid border -1\n";
    gp << "set boxwidth 0.5\n";
    gp << "plot '"<< aDataPath << "' "<< "using (100*$2/($2+$3+$4+$5+$6)):xtic(1) title columnheader(2) lc rgb \"blue\"\
      ,'' using (100*$3/($2+$3+$4+$5+$6)) title columnheader(3) lc rgb \"yellow\"\
        ,'' using (100*$4/($2+$3+$4+$5+$6)) title columnheader(4) lc rgb \"purple\"\
          ,'' using (100*$5/($2+$3+$4+$5+$6)) title columnheader(5) lc rgb \"orange\"\
            ,'' using (100*$6/($2+$3+$4+$5+$6)) title columnheader(6) lc rgb \"black\"";
  }
}

void CPlotter::plotTwoBarsGraph(const string& aDataPath, const string &aFilename, FormatSpec aFormat) {
  if ( !boost::filesystem::exists(aDataPath)) {
    cout << "Warning: no data to create plot: " << aFormat.title << endl;
  } else {
    Gnuplot gp;
    preSettings(gp, aFilename, aFormat, "bmargin");
    gp << "set style histogram cluster gap 1\n";
    gp << "set style data histograms\n";
    gp << "set style fill solid border -1\n";
    gp << "set yrange [0:*]\n";
    gp << "set ytics nomirror\n";
    gp << "set y2range [0:*]\n";
    gp << "set y2tics\n";
    gp << "plot '" << aDataPath << "' " << "u 2 axis x1y1 title columnheader(2) lc rgb \"blue\", '' u 3:xtic(1) axis x1y2 title columnheader(3) lc rgb \"yellow\"\n";
  }
}
void CPlotter::plotWeightedDistributionGraph(const string& aDataPath, const string &aFilename, FormatSpec aFormat) {
  if ( !boost::filesystem::exists(aDataPath)) {
    cout << "Warning: no data to create plot: " << aFormat.title << endl;
  } else {
    Gnuplot gp;
    preSettings(gp, aFilename, aFormat, "bmargin");
    gp << "set style data histograms\n";
    gp << "set style fill solid 1.00 border\n";
    gp << "set yrange [0:*]\n";
    gp << "set ytics nomirror\n";
    gp << "set y2range [0:*]\n";
    gp << "set y2tics\n";
    gp << "set xtic 4\n";
    gp << "set xtics rotate by -45\n";
    gp << "set boxwidth 0.005\n";
    gp << "everyfourth(col) = (int(column(col))%4 ==0)?stringcolumn(1):\"\"\n";
    gp << "plot '" << aDataPath << "' " << "using 2 axis x1y1 ti columnheader(2)   lc rgb \"blue\", '' u  3:xticlabels(everyfourth(1)) axis x1y2 with lines ti columnheader(3) lc rgb \"orange\"\n";
  }
}
