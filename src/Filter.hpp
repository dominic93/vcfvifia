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
 * @file Filter.hpp
 * VariantFilter wrapper
 */

#ifndef FILTER_HPP
#define FILTER_HPP

#include "Variant.h"

/**
 * enum for variant types
 */
enum EVariantType {
  SNP, /**< single-nucleotide polymorphism */
  MNP, /**< multi-nucleotide polymorphism */
  INS, /**< insertion */
  DEL, /**< deletion */
  COM, /**< complex (meaning any other) */
  DEFAULT /**< for computation purpose only */
};

/**
 * filter class
 */
class CFilter{
public:
  /**
   * strucht holding predefined filters
   */
  struct PredefinedFilters {
    int homozygosity = 0; /** 0 = no filter, 1 = homozygous variants only, -1 = heterozygous variants only **/
    int phased = 0; /** 0 = no filter, 1 = phased genotype only, -1 = unphased genotype only **/
    vector<EVariantType> variantTypeFilters; /** empty = no filter, otherwise: ignoring all variants with type not present in the vector **/
  };

  /**
   * constructor of the filter class
   *  @param aInfoFilter info filter string
   *  @param aInfoVariables info field variables
   *  @param aSampleFilter sample filter string
   *  @param aSampleVariables sample field variables
   *  @param aEffectFilter effect filter string
   *  @param aPredefinedFilters predefined filters
   */
  CFilter(std::string &aInfoFilter, std::map<std::string, vcflib::VariantFieldType> &aInfoVariables, std::string &aSampleFilter, std::map<std::string, vcflib::VariantFieldType> &aSampleVariables, std::string &aEffectFilter, PredefinedFilters aPredefinedFilters);

  /**
   * constructor of the filter class
   * info, sample and effect filter will always pass
   * @param aPredefinedFilters predefined filters
   */
  CFilter(PredefinedFilters aPredefinedFilters);

  /**
   * destructor of the filter class
   */
  virtual ~CFilter();

  /**
   * function to check if sample filter is passed
   * @param aVariant variant to check
   * @param aSample sample to check
   * @return true in case of passing, false otherwise
   */
  bool passesSampleFilter(vcflib::Variant &aVariant, std::string &aSample);

  /**
   * function to check if info filter is passed
   * @param aVariant variant to check
   * @return true in case of passing, false otherwise
   */
  bool passesInfoFilter(vcflib::Variant &aVariant);

  /**
   * function to check if effect filter is passed
   * @param aVariant variant to check
   * @return true in case of passing, false otherwise
   */
  bool passesEffectFilter(vcflib::Variant &aVariant);

  /**
   * function to check if predefined filters are passed
   * @param aVariant variant type of the record
   * @param aGenotypeSample genotype of first sample to check
   * @return true in case of passing, false otherwise
   */
  bool passesPredefinedFilters(EVariantType aVariantType, std::string &aGenotype);

private:
  /**
   * enum for filter behaviour
   */
  enum class EFilterBehaviour{USE, FAILING, PASSING};

  /**
   * enum used to determine sample filter behaviour
   * @see passesSampleFilter
   */
  EFilterBehaviour mSampleBehaviour;

  /**
   * enum used to determine info filter behaviour
   * @see passesInfoFilter
   */
  EFilterBehaviour mInfoBehaviour;

  /**
   * enum used to determine effect filter behaviour
   * @see passesEffectFilter
   */
  EFilterBehaviour mEffectBehaviour;

  /**
   * effect filter string
   * @see passesEffectFilter
   */
  std::string mEffectFilter;

  /**
   * predefined filters
   * @see PredefinedFilters
   */
  PredefinedFilters mPredefinedFilters;

  /**
   * pointer to info filter object
   * null in case of no filter specified
   * @see CFilter()
   */
 vcflib::VariantFilter *mpInfoFilter;

  /**
   * pointer to sample filter object
   * null in case of no filter specified
   * @see CFilter()
   */
 vcflib::VariantFilter *mpSampleFilter;

  /**
   * checks wheter the variables in the filter string are defined
   * @param aFilter the filter
   * @param aDefines the definitions from the header
   * @return true if filter is valid, false otherwise
   */
  bool checkFilterValidity(const std::string &aFilter, const std::map<std::string, vcflib::VariantFieldType> &aDefines);

  /**
   * checks wheter the genotype is homozygous
   * @param aGenotype string representation of genotype
   * @return true in case of homozygosity, false otherwise
   */
  bool isHomozygous(const std::string &aGenotype);

  /**
   * checks wheter the genotype is phased
   * @param aGenotype string representation of genotype
   * @return true in case of phased genotype, false otherwise
   */
  bool isPhased(const std::string &aGenotype);
};
#endif
