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

#include "Validator.hpp"
#include "utils.hpp"

#include <assert.h>
#include <iomanip>

void CValidator::computeValidationResults(const CBasicStatistic::CompareRes &aCompRes, EVariantType aVariantType, map<string, string> aGenotypesTruth, map<string, string> aGenotypesPrediction) {
  for(auto itr = aCompRes.both.begin(); itr != aCompRes.both.end(); ++itr) {
    const string &genotype = aGenotypesTruth[*itr];
    mValidation[aVariantType][genotype].gtCount++; //increment truth genotype
    mValidation[aVariantType][genotype].predictedCounts[aGenotypesPrediction[*itr]]++; //increment predicted genotype
  }
}

void CValidator::computeValidationResults(const CBasicStatistic::CompareDELRes &aCompDelRes) {
  assert(aCompDelRes.bothFromFirst.size() == aCompDelRes.bothFromSecond.size());
  for(unsigned int i = 0; i < aCompDelRes.bothFromFirst.size(); ++i) {
    string genotype = aCompDelRes.bothFromFirst.at(i).genotype;
    mValidation[DEL][genotype].gtCount++; //increment truth genotype
    mValidation[DEL][genotype].predictedCounts[aCompDelRes.bothFromSecond.at(i).genotype]++; //increment predicted genotype
  }
}

void CValidator::computeValidationResults(const CBasicStatistic::CompareINSRes &aCompInsRes) {
  assert(aCompInsRes.bothFromFirst.size() == aCompInsRes.bothFromSecond.size());
  for(unsigned int i = 0; i < aCompInsRes.bothFromFirst.size(); ++i) {
    string genotype = aCompInsRes.bothFromFirst.at(i).genotype;
    mValidation[INS][genotype].gtCount++; //increment truth genotype
    mValidation[INS][genotype].predictedCounts[aCompInsRes.bothFromSecond.at(i).genotype]++; //increment predicted genotype
  }
}

void CValidator::printValidationResults(const string &aFilenameTruth, const string &aFilenamePrediction) {
  string filenameTruth = aFilenameTruth.substr(aFilenameTruth.find_last_of("\\/")+1);
  string filenamePrediction = aFilenamePrediction.substr(aFilenamePrediction.find_last_of("\\/")+1);

  cout << "-- printing validation results --" << endl;
  for(auto itr = mValidation.begin(); itr != mValidation.end(); ++itr) {
    auto vec = predGenotypesForVariant(itr->first);
    if(vec.size() == 0) {
      continue;
    }
    cout << "-- " << CBasicStatistic::getNameForVariantType(itr->first) << " --" << endl;
    cout << filenameTruth << " / " << filenamePrediction << " " << utils::toString(vec) << endl;
    for(auto itr1 = itr->second.begin(); itr1 != itr->second.end(); ++itr1) {
      cout << itr1->first << " ";
      for(auto itr2 = vec.begin(); itr2 != vec.end(); ++itr2) {
        double total = itr1->second.gtCount;
        double pred = itr1->second.predictedCounts[*itr2];
        cout << fixed << setprecision(4) << pred/total << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
}

vector<string> CValidator::predGenotypesForVariant(EVariantType aVariantType) {
  vector<string> result;
  auto itr = mValidation[aVariantType].begin();
  for(; itr != mValidation[aVariantType].end(); ++itr) {
    auto itr2 = itr->second.predictedCounts.begin();
    for(; itr2 != itr->second.predictedCounts.end(); ++itr2) {
      result.push_back(itr2->first);
    }
  }
  sort( result.begin(), result.end() );
  result.erase( unique( result.begin(), result.end() ), result.end() );
  return result;
}
