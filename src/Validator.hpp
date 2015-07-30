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
 * @file Validator.hpp
 * functions to validate the results of a variant caller (kind of genotype concordance)
 */

#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP

#include "BasicStatistic.hpp"

/**
 * validator class
 */
class CValidator {
public:
  /**
   * function to compute validation results from comparison results
   * @param aCompRes compare results @see CBasicStatistic#CompareRes
   * @param aVariantType the type of the compare results @see EVariantType
   */ //TODO adjust comment
  void computeValidationResults(const CBasicStatistic::CompareRes &aCompRes, EVariantType aVariantType, std::map<std::string, std::string> aGenotypesTruth, std::map<std::string, std::string> aGenotypesPrediction);

  /**
   * function to compute validation results from deletion comparison results
   * @param aCompDelRes deletion compare results @see CBasicStatistic#CompareDELRes
   */ //TODO adjust comment
  void computeValidationResults(const CBasicStatistic::CompareDELRes &aCompDelRes);

  /**
   * function to compute validation results from insertion comparison results
   * @param aCompInsRes insertion compare results @see CBasicStatistic#CompareINSRes
   */ //TODO adjust comment
  void computeValidationResults(const CBasicStatistic::CompareINSRes &aCompInsRes);

  /**
   * function to print validation results
   * @param aFilenameTruth name of the reference file
   * @param aFilenamePrediction name of the prediction file
   */
  void printValidationResults(const std::string &aFilenameTruth, const std::string &aFilenamePrediction);

private:
  /**
   * struct holding validation results for a genotype
   * @see computeValidationResults
   */
  struct GtValidationRes {
    std::map<std::string, unsigned long> predictedCounts; /**< mapping genotypes to their prediction count> */
    unsigned long gtCount; /**< total count of the genotype */
  };

  /**
   * mapping variant types to mapping of genotype to genotype validation results @see GtValidationRes
   */
  std::map<EVariantType, std::map<std::string, GtValidationRes> > mValidation;

  /**
   * function to get the genotypes from prediction (by variant)
   * @param aVariantType the variant type
   * @return vector with string representation of the genotypes
   */
  std::vector<std::string> predGenotypesForVariant(EVariantType aVariantType);
};
#endif //VALIDATOR_HPP
