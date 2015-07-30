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

#include "Filter.hpp"
#include "utils.hpp"

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace vcflib;

CFilter::CFilter(PredefinedFilters aPredefinedFilters) {
  mInfoBehaviour = EFilterBehaviour::PASSING;
  mpInfoFilter = NULL;
  mSampleBehaviour = EFilterBehaviour::PASSING;
  mpSampleFilter = NULL;
  mEffectBehaviour = EFilterBehaviour::PASSING;
  mPredefinedFilters = aPredefinedFilters;
}

CFilter::CFilter(string &aInfoFilter, map<string, VariantFieldType> &aInfoVariables, string &aSampleFilter, map<string, VariantFieldType> &aSampleVariables, string &aEffectFilter, PredefinedFilters aPredefinedFilters) {
  mPredefinedFilters = aPredefinedFilters;
  mEffectFilter = aEffectFilter;
  if(mEffectFilter == "") {
    mEffectBehaviour = EFilterBehaviour::PASSING;
  } else {
    mEffectBehaviour = EFilterBehaviour::USE;
  }
  //Determine filter behaviour of info filter
  if(aInfoFilter != "") {
    if(checkFilterValidity(aInfoFilter, aInfoVariables)) {
      mpInfoFilter = new VariantFilter(aInfoFilter, VariantFilter::RECORD, aInfoVariables);
      mInfoBehaviour = EFilterBehaviour::USE;
    } else {
      mInfoBehaviour = EFilterBehaviour::FAILING;
      mpInfoFilter = NULL;
    }
  } else {
    mInfoBehaviour = EFilterBehaviour::PASSING;
    mpInfoFilter = NULL;
  }

  //Determine filter behvaviour of sample filter
  if(aSampleFilter != "") {
    if(checkFilterValidity(aSampleFilter, aSampleVariables)) {
      mpSampleFilter = new VariantFilter(aSampleFilter, VariantFilter::SAMPLE, aSampleVariables);
      mSampleBehaviour = EFilterBehaviour::USE;
    } else {
      mSampleBehaviour = EFilterBehaviour::FAILING;
      mpSampleFilter = NULL;
    }
  } else {
    mSampleBehaviour = EFilterBehaviour::PASSING;
    mpSampleFilter = NULL;
  }
}

CFilter::~CFilter() {
  delete mpInfoFilter;
  mpInfoFilter = NULL;
  delete mpSampleFilter;
  mpSampleFilter = NULL;
}

bool CFilter::passesPredefinedFilters(EVariantType aVariantType, string &aGenotype) {
  if(mPredefinedFilters.phased != 0 || mPredefinedFilters.homozygosity != 0) {
    bool homozygous = isHomozygous(aGenotype);
    bool phased = isPhased(aGenotype);
    if((homozygous && mPredefinedFilters.homozygosity == -1) ||
       (!homozygous && mPredefinedFilters.homozygosity == 1)) {
      return false;
    }
   if((phased && mPredefinedFilters.phased == -1) ||
       (!phased && mPredefinedFilters.phased == 1)) {
      return false;
    }
  }
  auto typeFilters = mPredefinedFilters.variantTypeFilters;
  bool result = true;
  if(!typeFilters.empty()) {
    result = false;
    for(auto itr = typeFilters.begin(); itr != typeFilters.end(); ++itr) {
      if((*itr) == aVariantType) {
        return true;
      }
    }
  }
  return result;
}

bool CFilter::passesSampleFilter(Variant &aVariant, string &aSample) {
  switch(mSampleBehaviour) {
    case EFilterBehaviour::PASSING:
      return true;
      break;
    case EFilterBehaviour::FAILING:
      return false;
      break;
    case EFilterBehaviour::USE:
      return mpSampleFilter->passes(aVariant, aSample);
      break;
    default:
      return false;
    }
}

bool CFilter::passesEffectFilter(Variant &aVariant) {
  if(mEffectBehaviour == EFilterBehaviour::PASSING) {
    return true;
  }
  for(auto itr = aVariant.info.begin(); itr != aVariant.info.end(); ++itr) {
    for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {
      size_t found = (*itr2).find(mEffectFilter);
      if (found!=std::string::npos) {
        return true;
      }
    }
  }
  return false;
}

bool CFilter::passesInfoFilter(Variant &aVariant) {
  string sample = "";
  switch(mInfoBehaviour) {
    case EFilterBehaviour::PASSING:
      return true;
      break;
    case EFilterBehaviour::FAILING:
      return false;
      break;
    case EFilterBehaviour::USE:
      return mpInfoFilter->passes(aVariant, sample);
      break;
    default:
      return false;
  }
}

bool CFilter::checkFilterValidity(const string &aFilter, const map<string, VariantFieldType> &aDefines) {
  bool result = true;
  string cleanedFilter = aFilter;
  vector<string> keyWords = {"(", ")", "&", "|", "<", ">", "=", "!"};
  for(auto itr = keyWords.begin(); itr != keyWords.end(); ++itr) {
    boost::erase_all(cleanedFilter, *itr);
  }
  vector<string> tokens;
  vector<string> cleanedTokens;
  boost::split(tokens, cleanedFilter, boost::is_any_of(" "), boost::token_compress_on);
  for(auto tokenItr = tokens.begin(); tokenItr != tokens.end(); ++tokenItr) {
    if(!utils::firstIsDigit(*tokenItr) && !(*tokenItr).empty()) {
      cleanedTokens.push_back(*tokenItr);
    }
  }

  for(auto itr = cleanedTokens.begin(); itr != cleanedTokens.end(); ++itr) {
    bool found = false;
    for(auto itr1 = aDefines.begin(); itr1 != aDefines.end(); ++itr1) {
      if((*itr) == (itr1->first)) {
        found = true;
        break;
      }
    }
    if(!found) {
      cerr << "Warning: couldn't find specification of " << (*itr) << " in header of vcf file" << endl;
      cerr << "Default behaviour: filter will always fail!" << endl;
      result = false;
    }
  }
  return result;
}

bool CFilter::isHomozygous(const string &aGenotype) {
  vector<string> genotype = split(aGenotype, "|/");
  if(atoi(genotype[0].c_str()) == atoi(genotype[1].c_str())) {
    if(genotype[0] != ".") {
      return true;
    }
  }
  return false;
}

bool CFilter::isPhased(const string &aGenotype) {
  vector<string> gtPhased = split(aGenotype, "|");
  if(gtPhased.size() == 2) {
    return true;
  }
  return false;
}
