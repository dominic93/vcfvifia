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
* @file main.cpp
* main logic of VCFvifia (including argument parsing)
*/

#include "Plotter.hpp"
#include "Variant.h"
#include "BasicStatistic.hpp"
#include "Filter.hpp"
#include "Validator.hpp"
#include "config.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <exception>
#include <vector>
#include <string>

using namespace std;
namespace po = boost::program_options;
namespace bf = boost::filesystem;

int main( int argc, const char* argv[] ) {

  //Default options
  string plotFormat = "pdf";
  bool writeCompResToFiles = false;
  bool plot = true;
  bool countGenotypeOnly = false;
  bool validate = false;
  bool printOut = false;

  if(!isGnuplotInstalled()) {
    cout << "Warning: Plots are disabled because gnuplot is not found." << endl;
    cout << "Ensure gnuplot is installed on your system to use the plotting feature." << endl;
    cout << endl;
    plot = false;
  }

  //Describe allowed options
  CFilter::PredefinedFilters predefFilters;
  CBasicStatistic::Thresholds thresholds;
  vector<string> vcfFilenames;
  vector<string> variantFilters;
  string infoFilter, infoFilter2;
  string effectFilter, effectFilter2;
  string sampleFilter, sampleFilter2;
  po::options_description generalOptions("General options");
  generalOptions.add_options()
    ("vcfFiles,v", po::value< vector<string> >(&vcfFilenames)->multitoken()->required(), "VCF files to be processed (at least one file required).\n")
    ("help,h", "Print help.");

  po::options_description plotOptions("Plotting options");
  plotOptions.add_options()
    ("png", "Create plots in png format (default: pdf).\n")
    ("noplot","Disable plots.");

  po::options_description filterOptions("Filter options");
  filterOptions.add_options()
    ("infoFilter", po::value< string >(&infoFilter), "Do statistics using this info filter for the first vcf file (see gtFilter for further information).\n")
    ("gtFilter", po::value< string >(&sampleFilter), "Do statistics using this sample filter for the first vcf file.\nAllowed operators: =, !, <, >, |, &\nEnsure that every logical unit is divided by a space (e.g. : \"! ( G = 1|1 )\", but not: \"!(GT = 1|1)\".\n")
    ("effFilter", po::value< string >(&effectFilter), "Do statistics using this effect filter for the first vcf file (use this filter if info filter can't find annotation, e.g. annotations from SnpEff or VEP)\n")
    ("infoFilter2", po::value< string >(&infoFilter2), "Info filter for the second vcf file (see gtFilter for further information).\n")
    ("gtFilter2", po::value< string >(&sampleFilter2), "Sample filter for the second vcf file (see gtFilter for further information).\n")
    ("effFilter2", po::value< string >(&effectFilter2), "Do statistics using this effect filter for the first vcf file (use this filter if info filter can't find annotation, e.g. annotations from SnpEff or VEP)\n")
    ("variants", po::value< vector<string> >(&variantFilters)->multitoken(), "Do statistics regarding only variants of specified types.\nVariant types can be: snp, mnp, insertion, deletion, other. Default: all\n")
    ("phasedOnly", "Records with unphased genotypes will be filtered out.\n")
    ("unphasedOnly", "Records with phased geontypes will be filtered out.\n")
    ("homozygousOnly", "Records with heterozygous genotypes will be filtered out.\n")
    ("heterozygousOnly", "Records with homozygous genotypes will be filtered out.\n")
    ("output", "Print records which passed all filters (or each record in case of no filters).\n")
    ("gtOnly", "Print genotype counts only.");

  po::options_description compareOptions("Comparison options");
  compareOptions.add_options()
    ("comp", "Writing detailed comparison information to files.\n")
    ("delLenDiff", po::value< unsigned int >(&thresholds.maxDelLengthDiff), "Maximum allowed length difference of two deletions to be detected to be identical.\n")
    ("delCenDiff", po::value< unsigned int >(&thresholds.maxDelCenterDiff), "Maximum allowed difference of center points of two deletions to be detected to be idential.\n")
    ("insLenDiff", po::value< unsigned int >(&thresholds.maxInsLengthDiff), "Maximum allowed length difference of two insertions to be detected to be identical.\n")
    ("insPosDiff", po::value< unsigned int >(&thresholds.maxInsPosDiff), "Maximum allowed difference of the positions of two insertions to be detected to be idential.\n")
    ("validate", "Printing genotype concordance results, where the first file is regarded as the reference and the second as the prediction.");

  po::options_description allOptions("Allowed options");
  if(plot) {
    allOptions.add(generalOptions).add(plotOptions).add(filterOptions).add(compareOptions);
  } else {
    allOptions.add(generalOptions).add(filterOptions).add(compareOptions);
  }

  //Parse options
  po::variables_map options;
    try {
      po::store(po::parse_command_line(argc, argv, allOptions), options);
      po::notify(options);
    } catch(std::exception& e) {
      cerr << "Error: " << e.what() << endl;
      cerr << allOptions << endl;
      return 1;
    }

  //Set flags depending on parsed options
  if(options.count("output")) {
    printOut = true;
  }
  if(options.count("validate")) {
    validate = true;
  }
  if(options.count("gtOnly")) {
    countGenotypeOnly = true;
  }
  if(options.count("comp")) {
    writeCompResToFiles = true;
  }
  if(options.count("noplot")) {
    plot = false;
  }
  if(options.count("help")) {
    cerr << allOptions << endl;
    return 1;
  }
  if(options.count("png")) {
    plotFormat = "png";
  }
  if(options.count("phasedOnly")) {
    predefFilters.phased = 1;
  }
  if(options.count("unphasedOnly")) {
    predefFilters.phased = -1;
  }
  if(options.count("homozygousOnly")) {
    predefFilters.homozygosity = 1;
  }
  if(options.count("heterozygousOnly")) {
    predefFilters.homozygosity = -1;
  }

  if(!variantFilters.empty()) {
    for(auto varFilter = variantFilters.begin(); varFilter != variantFilters.end(); ++varFilter) {
      if(variantTypes.count(*varFilter)) {
        predefFilters.variantTypeFilters.push_back(variantTypes[*varFilter]);
      }
    }
  }

  //Parse vcf files
  unsigned int vcfFileCount = vcfFilenames.size();
  vcflib::VariantCallFile *vcfFiles = new vcflib::VariantCallFile[vcfFileCount];
  for(unsigned int i = 0; i < vcfFileCount; ++i) {
    string currentFile = vcfFilenames[i];
    try{
      vcfFiles[i].open(currentFile);
    } catch(std::exception& e) {
      cerr << "Error: cannot open vcf_file: " << currentFile << endl;
      cerr << "description: "<< e.what() << endl;
      cerr << allOptions << endl;
      delete[]vcfFiles;
      return 1;
    } catch(...){
      cerr << "Error: cannot open vcf_file: " << currentFile << endl;
      cerr << allOptions << endl;
      delete[]vcfFiles;
      return 1;
    }

    if(!vcfFiles[i].is_open()){
      cerr << "Error: cannot open vcf_file: " << currentFile << endl;
      delete[]vcfFiles;
      return 1;
    }
  }

  cout << endl;
  cout << "--   printing basic statistics    --" << endl;
  cout << endl;

  CPlotter plotter(plotFormat);
  CBasicStatistic basicStatistic(thresholds);

  //main logic
  for(unsigned int i = 0; i < vcfFileCount; ++i) {
    cout << "-- file " << vcfFilenames[i] << endl;
    CFilter *filter;
    if(i == 0) {
      filter = new CFilter(infoFilter, vcfFiles[i].infoTypes, sampleFilter, vcfFiles[i].formatTypes,effectFilter, predefFilters);
    } else if (i == 1){
      filter = new CFilter(infoFilter2, vcfFiles[i].infoTypes, sampleFilter2, vcfFiles[i].formatTypes, effectFilter2, predefFilters);
    } else {
      filter = new CFilter(predefFilters);
    }
    basicStatistic.doBasicStats(vcfFiles[i], vcfFilenames[i], filter, printOut);
    if(countGenotypeOnly) {
      basicStatistic.printGenotype(vcfFilenames[i]);
    } else {
      basicStatistic.printBasicStats(vcfFilenames[i]);
      if(plot) {
        CPlotter::FormatSpec formatSubst("Substitution types", "Substitution", "Counts", "", "substTypes");
        CPlotter::FormatSpec formatIndel("Indel length distribution", "Length", "Counts", "", "indelLengthDist");
        CPlotter::FormatSpec formatStats("Basic statistics", "", "\% of total", "", "basicStats");
        CPlotter::FormatSpec formatSubstFreq("Allele frequency weighted substitution types", "Substitution","Counts", "Allele frequency",  "AFweightedSubstTypes");
        CPlotter::FormatSpec formatIndelFreq("Allele frequency weighted Indel length distribution", "Length", "Counts", "Allele frequency", "AFweightedIndelLengthDist");
        bf::path tempDir = bf::unique_path();
        bf::path tempFile(tempDir /"temp.dat");
        bf::path tempFile2(tempDir /"temp2.dat");
        bf::path tempFile3(tempDir /"temp3.dat");
        if(!bf::create_directory(tempDir)) {
          cerr << "Warning: cannot create basic statistics for file: " << vcfFilenames[i] << endl;
          cerr << "Could not create temporary directory '" << tempDir.c_str() << "' needed to print basic statistics" << endl;
          delete filter;
          continue;
        }
        basicStatistic.writeBasicStatsToFile(vcfFilenames[i], tempFile.c_str());
        basicStatistic.writeSubstFreqDataToFile(vcfFilenames[i], tempFile2.c_str());
        basicStatistic.writeIndelFreqDataToFile(vcfFilenames[i], tempFile3.c_str());
        plotter.plotStackedBarGraph(tempFile.c_str(), vcfFilenames[i], formatStats);
        plotter.plotBarGraph(basicStatistic.getSubstData(vcfFilenames[i]), vcfFilenames[i], formatSubst);
        plotter.plotTwoBarsGraph(tempFile2.c_str(), vcfFilenames[i], formatSubstFreq);
        plotter.plotDistributionGraph(basicStatistic.getIndelData(vcfFilenames[i]), vcfFilenames[i], formatIndel);
        plotter.plotWeightedDistributionGraph(tempFile3.c_str(), vcfFilenames[i], formatIndelFreq);
        remove_all(tempDir); //Clean up
        cout << endl;
      }
    }
    delete filter;
  }
  cout << endl;

  if(vcfFileCount > 1) { //Compare files (if there are at least two)
    if(validate) {
      CValidator validator;
      basicStatistic.compareTwo(vcfFilenames[0], vcfFilenames[1], writeCompResToFiles, validate);
      vector<EVariantType> variants = {SNP, MNP, COM};
      for(auto itr = variants.begin(); itr != variants.end(); ++itr) {
        validator.computeValidationResults(basicStatistic.getCompareResByType(*itr), *itr, basicStatistic.getVreprToGT(vcfFilenames[0]), basicStatistic.getVreprToGT(vcfFilenames[1]));
      }
      validator.computeValidationResults(basicStatistic.getCompareINSRes());
      validator.computeValidationResults(basicStatistic.getCompareDELRes());

      validator.printValidationResults(vcfFilenames[0], vcfFilenames[1]);
    } else {
      basicStatistic.compareTwo(vcfFilenames[0], vcfFilenames[1], writeCompResToFiles, validate);
    }
  }
  delete[]vcfFiles;
  return 0;
}
