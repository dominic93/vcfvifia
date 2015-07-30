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
 * @file BasicStatistic.hpp
 * functions to perform basic statistic and access the produced data
 */

#ifndef BASICSTATISTIC_HPP
#define BASICSTATISTIC_HPP

#include "Filter.hpp"
#include "utils.hpp"

#include "Variant.h"
#include "split.h"
#include "convert.h"

#include <getopt.h>
#include <boost/filesystem.hpp>
#include <fstream>

/**
 * substitution types
 */
enum ESubstitutionType {
  AT,
  AG,
  AC,
  TA,
  TG,
  TC,
  GA,
  GT,
  GC,
  CA,
  CG,
  CT
};

/**
 * string representations of the variant types
 * @see EVariantType
 */
static const char* variantNames[DEFAULT+1] = { "SNP", "MNP", "INSERTION", "DELETION", "COMPLEX", "NONE" };

/**
 * mapping strings to variant types
 */
static std::map <std::string, EVariantType> variantTypes = {{"snp", SNP}, {"mnp", MNP}, {"insertion", INS}, {"deletion", DEL}, {"other", COM}};

/**
 * string representations of the substitution types
 * @see ESubstitutionType
 */
static const char* substNames[CT+1] = { "A>T", "A>G", "A>C", "T>A", "T>G", "T>C", "G>A", "G>T", "G>C", "C>A", "C>G", "C>T" };

/**
 * mapping strings to substitution types
 */
static std::map <std::string, ESubstitutionType> substitutionTypes = {{"AT", AT},
                                                                  {"AG", AG},
                                                                  {"AC", AC},
                                                                  {"TA", TA},
                                                                  {"TG", TG},
                                                                  {"TC", TC},
                                                                  {"GA", GA},
                                                                  {"GT", GT},
                                                                  {"GC", GC},
                                                                  {"CA", CA},
                                                                  {"CG", CG},
                                                                  {"CT", CT}};

/**
 * basic statistic class
 */
class CBasicStatistic {
public:
  /**
   * thresholds for the comparison of insertions and deletions
   */
  struct Thresholds {
    unsigned int maxDelLengthDiff = 0;
    unsigned int maxDelCenterDiff = 0;
    unsigned int maxInsLengthDiff = 0;
    unsigned int maxInsPosDiff = 0;
  };

  /**
   * constructor of the basic statistic class
   */
  CBasicStatistic(Thresholds aThresholds);

  /**
   * comparison results
   * @see computeComparison @see CValidator#computeValidationResults
   */
  struct CompareRes {
    std::vector<std::string> firstOnly;
    std::vector<std::string> secondOnly;
    std::vector<std::string> both;
  };

  /**
   * representation of insertion for less restrictive comparison
   */
  struct InsRep {
    std::string chrom; /**< the chromosome */
    std::string genotype; /**< the genotype */
    int pos; /**< the position on the chromosome */
    int length; /**< the length of the insertion */
    string str() const{
      stringstream sstream;
      sstream << chrom;
      sstream << " ";
      sstream << pos;
      sstream << " ";
      sstream << length;
      sstream << " ";
      sstream << genotype;
      return sstream.str();
    }
  };

  /**
   * representation of deletion for less restrictive comparison
   */
  struct DelRep {
    std::string chrom; /**< the chromosome */
    std::string genotype; /**< the genotype */
    int center; /**< the center of the deletion */
    int pos; /**< the position on the chromosome (for result printing only)*/
    int length; /**< the length of the deletion */
    string str() const {
      stringstream sstream;
      sstream << chrom;
      sstream << " ";
      sstream << pos;
      sstream << " ";
      sstream << length;
      sstream << " ";
      sstream << genotype;
      return sstream.str();
    }
  };

  /**
   * insertion comparison results
   * @see computeComparison @see CValidator#computeValidationResults
   */
  struct CompareINSRes {
    std::vector<InsRep> firstOnly;
    std::vector<InsRep> secondOnly;
    std::vector<InsRep> bothFromFirst;
    std::vector<InsRep> bothFromSecond;
  };

  /**
   * deletion comparison results
   * @see computeComparison @see CValidator#computeValidationResults
   */
  struct CompareDELRes {
    std::vector<DelRep> firstOnly;
    std::vector<DelRep> secondOnly;
    std::vector<DelRep> bothFromFirst;
    std::vector<DelRep> bothFromSecond;
  };


  /*
   * getter for mapping from vrepr to gt
   * @param aFilename
   * @return mapping from vrepr to gt
   */
   std::map<std::string, std::string> getVreprToGT(const std::string &aFilename);

  /**
   * function to do basic statistics,
   * counting variants, assigning a variant type to each record, collect information used by other functions (= should always be called at first)
   * @param aVariantFile vcf file to be processed
   * @param aFilename name of the vcf file (used internally for storing collected data)
   * @param apFilter filter object
   * @param aPrintOut if true, each record who passed all filters will be printed out
   */
  void doBasicStats(vcflib::VariantCallFile &aVariantFile, const std::string &aFilename, CFilter *apFilter, bool aPrintOut);

  /**
   * function to access comparison results by type @see CompareRes
   * @param aVariantType the variant type @see EVariantType
   * @return CompareRes for type
   */
  inline CompareRes getCompareResByType(EVariantType aVariantType) {
    return mCompareResByType[aVariantType];
  }

  /**
   * function to access insertion comparison results @see CompareINSRes
   * @return CompareINSRes
   */
  inline CompareINSRes getCompareINSRes() {
    return mCompareINSRes;
  }

  /**
   * function to access deletion comparison results @see CompareDELRes
   * @return CompareDELRes
   */
  inline CompareDELRes getCompareDELRes() {
    return mCompareDELRes;
  }

  /**
   * function to compare two vcf files
   * computes the intersection and the complements
   * assumes that the files have been processed using doBasicStats before
   * @param aFilenameOne name of the first file
   * @param aFilenameTwo name of the second file
   * @param aEval in case of true deeper statistic will be printed
   * @aCompToFile in case of true the sets will be saved in files
   */
  void compareTwo(const std::string &aFilenameOne, const std::string &aFilenameTwo, bool aCompToFile, bool aEval);

  /**
   * function to print genotype information
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   */
  void printGenotype(const std::string &aFilename);

  /**
   * function to print basic statistic
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   */
  void printBasicStats(const std::string &aFilename);

  /**
   * function to access substitution data
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   * @return vector with pairs of substitution type and count
   */
  std::vector<pair<std::string, unsigned int>>  getSubstData(const std::string &aFilename);

  /**
   * function to access indel distriution
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   * @return vector with pairs of indel length and count
   */
  std::vector<pair<int, unsigned int>> getIndelData(const std::string &aFilename);

  /**
   * function to write basic statistic to file
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   * @param aOutfile name of the output file
   */
  void writeBasicStatsToFile(const std::string &aFilename, const std::string &aOutfile);

  /**
   * function to write substitution with frequency data to file
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   * @return aOutfile name of the output file
   */
  void writeSubstFreqDataToFile(const std::string &aFilename, const std::string &aOutfile);

  /**
   * function to write indels with frequency data to file
   * assumes that the file has been processed using doBasicStats before
   * @param aFilename name of the vcf file
   * @return aOutfile name of the output file
   */
  void writeIndelFreqDataToFile(const std::string &aFilename, const std::string &aOutfile);

 /**
  * function to get the name for a variant type @see CFilter#EVariantType
  */
  static const char * getNameForVariantType(EVariantType aType) {
    if(aType < DEFAULT) {
      return variantNames[aType];
    } else {
      return "ERROR";
    }
  }

private:
  /**
   * struct holding basic statistics
   * @see doBasicStats
   */
  struct BasicStats {
    unsigned int passedFilter = 0;
    unsigned int records = 0;
    unsigned int samples = 0;
    unsigned int snps = 0;
    unsigned int mnps = 0;
    unsigned int ins = 0;
    unsigned int dels = 0;
    unsigned int largerIns = 0;
    unsigned int largerDels = 0;
    unsigned int complexVar = 0;
    unsigned int multiAllSNPs = 0;
    unsigned int multiAll = 0;
    unsigned int substTypes[CT+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double substFrequency[CT+1] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
    std::vector<std::string> snpsVec;
    std::vector<std::string> mnpsVec;
    std::vector<std::string> compVec;
    std::vector<InsRep> insRepVec;
    std::vector<DelRep> delRepVec;
    std::map<int, unsigned int> indelLengthCounts;
    std::map<int, double> indelLengthFrequency;
    std::map<std::string, unsigned int> genotypes;
    std::map<string, string> VrepToGT;
  };

  /**
  * the thresholds used in comparison of two InsRep/DelRep objects
  */
  Thresholds mThresholds;

  /**
   * mapping sequence name to int
   * this structure is used to realise a strict order of InsRep and DelRep objects
   */
  map<string, unsigned int> mOrderMap;

  /**
   * mapping variant types to compare res
   */
  std::map<EVariantType, CompareRes> mCompareResByType;

  /**
   * @see CompareINSRes
   */
  CompareINSRes mCompareINSRes;

  /**
   * @see CompareDELRes
   */
  CompareDELRes mCompareDELRes;

  /**
   * mapping filenames to BasicStats @see BasicStats
   */
  std::map<std::string, BasicStats> mStatsPerFile;

  /**
   * function to compute a records insertion representation
   * @param aVariant the record
   * @param aGenotype the genotype of the first sample
   * @return insertion representation @see InsRep
   */
  InsRep computeInsRep(vcflib::Variant aVariant, const string &aGenotype);

  /**
   * function to compute a records deletion representation
   * @param aVariant the record
   * @param aGenotype the genotype of the first sample
   * @return deletion representation @see DelRep
   */
  DelRep computeDelRep(vcflib::Variant aVariant, const string &aGenotype);

  /**
   * function to compute the intersection and complements of two vectors containing strings
   * @param aFirst first vector
   * @param aSecond second vector
   * @param aResult vectors the result will be stored in @see CompareRes
   */
  void computeComparison(std::vector<std::string> &aFirst, std::vector<std::string> &aSecond, CompareRes &aResult);

  /**
   * function to compute the intersection and complements of two vectors containing InsRep @see InsRep
   * @param aFirst first vector
   * @param aSecond second vector
   * @param aResult vectors the result will be stored in @see CompareRes
   */
  void computeInsComparison(std::vector<InsRep> &aFirst, std::vector<InsRep> aSecond, CompareRes &aResult);

  /**
   * function to compute the intersection and complements of two vectors containing DelRep @see DelRep
   * @param aFirst first vector
   * @param aSecond second vector
   * @param aResult vectors the result will be stored in @see CompareRes
   */
  void computeDelComparison(std::vector<DelRep> &aFirst, std::vector<DelRep> aSecond, CompareRes &aResult);

  /**
   * function to compare two vrepr strings without genotype
   * @param aFirst first vrepr
   * @param aSecond second vrepr
   * @return true in case of aFirst appears before aSecond in strict order, false otherwise
   */
  bool compareStringsWithoutGT(const std::string &aFirst, const std::string &aSecond);

  /**
   * function to compare two InsRep objects @see InsRep
   * @param aFirst first InsRep
   * @param aSecond second InsRep
   * @return true in case of aFirst appears before aSecond in strict order, false otherwise
   */
  bool compareInsRep(const InsRep &aFirst, const InsRep &aSecond);

  /**
   * function to compare two DelRep objects @see DelRep
   * @param aFirst first DelRep
   * @param aSecond second DelRep
   * @return true in case of aFirst appears before aSecond in strict order, false otherwise
   */
  bool compareDelRep(const DelRep &aFirst, const DelRep &aSecond);

  /**
   * function to set the larger number last in unphased genotype (1/0 -> 0/1)
   * phased genotypes will be ignored
   * @param aGenotype genotype string
   */
  void normalizeGenotype(std::string &aGenotype);

  /**
   * function to get the name for a substitution type @see ESubstitutionType
   */
  inline const char * getNameForSubstType(ESubstitutionType aType) {
    return substNames[aType];
  }

  /**
   * function to compute the allele frequency of a specified allele
   * @param aSamples the samples (where the genotype info is used from)
   * @param aAllIndex the index of the allele
   * @return allele frequency
   */
  double computeAlleleFreq(std::map<std::string, std::map<std::string, std::vector<std::string> > > aSamples, int aAllIndex);

  /**
   * function to compute and update the indel length + frequency information of a record
   * @param aRecord the record
   * @param aBasicStats the statistics to be updated
   */
  void updateIndelLengthFrequency(vcflib::Variant &aRecord, BasicStats &aBasicStats);

  /**
   * function to get the substition type
   * @param aFrom substitution from
   * @param aTo substitution to
   * @return @see ESubstitutionType the computed substitution type
   */
  inline ESubstitutionType getSubstitutionType(const std::string &aFrom, const std::string &aTo) {
    std::string substType = aFrom + aTo;
    return substitutionTypes[substType];
  }

  /**
   * function to get the variant type
   * @param aRef reference
   * @param aAlt alternative
   * @return @see EVariantType the computed variant type
   */
  inline  EVariantType getVariantType(const std::string &aRef, const std::string &aAlt) {
    if(isSNP(aRef, aAlt)) {
      return SNP;
    }
    if(isMNP(aRef, aAlt)) {
      return MNP;
    }
    if(isDEL(aRef, aAlt)) {
      return DEL;
    }
    if(isINS(aRef, aAlt)) {
      return INS;
    }
    return COM;
  }

  /**
   * function to check if variant is an insertion
   * @param aRef reference
   * @param aAlt alternative
   * @return true in case of an insertion, false otherwise
   */
  inline bool isINS(const std::string &aRef, const std::string &aAlt) {
    if (aRef.size() == 1  && aRef.size() < aAlt.size()) {
      if(aRef.compare(aAlt.substr(0, aRef.size())) == 0) {
        return true;
      }
      return false;
    }
    return false;
  }

  /**
   * function to check if variant is a deletion
   * @param aRef reference
   * @param aAlt alternative
   * @return true in case of a deletion, false otherwise
   */
  inline bool isDEL(const std::string &aRef, const std::string &aAlt) {
    if (aAlt.size() == 1 && aRef.size() > aAlt.size()) {
      if (aAlt.compare(aRef.substr(0, aAlt.size())) == 0) {
        return true;
      }
      return false;
    }
    return false;
  }

  /**
   * function to check if variant is a mnp
   * @param aRef reference
   * @param aAlt alternative
   * @return true in case of a mnp, false otherwise
   */
  inline bool isMNP(const std::string &aRef, const std::string &aAlt) {
    if (aRef.size() == aAlt.size()) {
      if (aRef.size() != 1) {
        return true;
      }
    }
    return false;
  }

  /**
   * function to check if variant is a snp
   * @param aRef reference
   * @param aAlt alternative
   * @return true in case of a snp, false otherwise
   */
  inline bool isSNP(const std::string &aRef, const std::string &aAlt) {
    if (aRef.size() == aAlt.size()) {
      if (aRef.size() == 1) {
        if((aRef[0] != '.') && (aAlt[0] != '.')) {
          return true;
        }
      }
    }
    return false;
  }
};

#endif //BASICSTATISTIC_HPP
