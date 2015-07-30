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

#include "BasicStatistic.hpp"

#include <iomanip>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace vcflib;

CBasicStatistic::CBasicStatistic(Thresholds aThresholds):mThresholds(aThresholds) { }

void CBasicStatistic::compareTwo(const string &aFilenameOne, const string &aFilenameTwo, bool aCompToFile, bool aEval) {
  BasicStats statsFile1 = mStatsPerFile[aFilenameOne];
  BasicStats statsFile2 = mStatsPerFile[aFilenameTwo];

  CompareRes &snps  = mCompareResByType[SNP];
  CompareRes &mnps  = mCompareResByType[MNP];
  CompareRes &ins   = mCompareResByType[INS];
  CompareRes &dels  = mCompareResByType[DEL];
  CompareRes &comps = mCompareResByType[COM];

  //Compute comparison
  computeComparison(statsFile1.snpsVec, statsFile2.snpsVec, snps);
  computeComparison(statsFile1.mnpsVec, statsFile2.mnpsVec, mnps);
  computeInsComparison(statsFile1.insRepVec, statsFile2.insRepVec, ins);
  computeDelComparison(statsFile1.delRepVec, statsFile2.delRepVec, dels);
  computeComparison(statsFile1.compVec, statsFile2.compVec, comps);

  //Print result summary
  string filename = aFilenameOne.substr(aFilenameOne.find_last_of("\\/")+1);
  string filename2 = aFilenameTwo.substr(aFilenameTwo.find_last_of("\\/")+1);

  cout << "-- printing comparison results --" << endl;
  cout << endl;

  cout << "-- file " << filename << " only --" << endl;
  cout << endl;
  cout << "# snps: " << snps.firstOnly.size() << endl;
  cout << "# mnps: " << mnps.firstOnly.size() << endl;
  cout << "# insertions: " << ins.firstOnly.size() << endl;
  cout << "# deletions: " << dels.firstOnly.size() << endl;
  cout << "# complex: " << comps.firstOnly.size() << endl;
  cout << endl;

  cout << "-- file " << filename2 << " only --" << endl;
  cout << endl;
  cout << "# snps: " << snps.secondOnly.size() << endl;
  cout << "# mnps: " << mnps.secondOnly.size() << endl;
  cout << "# insertions: " << ins.secondOnly.size() << endl;
  cout << "# deletions: " << dels.secondOnly.size() << endl;
  cout << "# complex: " << comps.secondOnly.size() << endl;
  cout << endl;

  cout << "-- file " << filename << " and " << filename2 << " --" << endl;
  cout << endl;
  cout << "# snps: " << snps.both.size() << endl;
  cout << "# mnps: " << mnps.both.size() << endl;
  cout << "# insertions: " << ins.both.size() << endl;
  cout << "# deletions: " << dels.both.size() << endl;
  cout << "# complex: " << comps.both.size() << endl;
  cout << endl;

  cout << "fraction of records in " << filename2 << " compared to " << filename << ": ";
  cout <<  fixed << setprecision(4) << ((double)statsFile2.passedFilter)/((double)statsFile1.passedFilter) << endl;
  cout << endl;

  if(aEval) {
    double truePositive = snps.both.size() + mnps.both.size() + ins.both.size() + dels.both.size() + comps.both.size();
    double falseNegative = snps.firstOnly.size() + mnps.firstOnly.size() + ins.firstOnly.size() + dels.firstOnly.size() + comps.firstOnly.size();
    double falsePositive = snps.secondOnly.size() + mnps.secondOnly.size() + ins.secondOnly.size() + dels.secondOnly.size() + comps.secondOnly.size();
    double recall = truePositive/(truePositive + falseNegative);
    double precision = truePositive/(truePositive + falsePositive);
  cout << "-- file " << filename << " is handled as truth and " << filename2 << " as prediction" << endl;
  cout << endl;
  cout << "recall:    " << fixed << setprecision(4) << recall << endl;
  cout << "precision: " << fixed << setprecision(4) << precision << endl;
  cout << endl;
  }

  if(aCompToFile) { //Write sets to file

    string filenameBoth = filename + "_and_" + filename2;

    string headerEndingOnly = " only";


    string filenameEnding = "_OnlySNPs.txt";

    string outfile, header;

    //SNPs
    outfile = filename + filenameEnding;
    header = "# SNPs in ";
    header += filename;
    header += " only";
    utils::printVectorToFile(snps.firstOnly, outfile, header);

    outfile = filename2 + filenameEnding;
    header = "# SNPs in ";
    header += filename2;
    header += " only";
    utils::printVectorToFile(snps.secondOnly, outfile, header);

    outfile = filenameBoth + "_SNPs.txt";
    header = "# SNPs in ";
    header += filenameBoth;
    utils::printVectorToFile(snps.both, outfile, header);

    //MNPs
    filenameEnding = "_OnlyMNPs.txt";
    outfile = filename + filenameEnding;
    header = "# MNPs in ";
    header += filename;
    header += " only";
    utils::printVectorToFile(mnps.firstOnly, outfile, header);

    outfile = filename2 + filenameEnding;
    header = "# MNPs in ";
    header += filename2;
    header += " only";
    utils:: printVectorToFile(mnps.secondOnly, outfile, header);

    outfile = filenameBoth + "_MNPs.txt";
    header = "# MNPs in ";
    header += filenameBoth;
    utils::printVectorToFile(mnps.both, outfile, header);

    //INSERTIONS
    filenameEnding = "_OnlyInserts.txt";
    outfile = filename + filenameEnding;
    header = "# Insertions in ";
    header += filename;
    header += " only";
    utils::printVectorToFile(ins.firstOnly, outfile, header);

    outfile = filename2 + filenameEnding;
    header = "# Insertions in ";
    header += filename2;
    header += " only";
    utils::printVectorToFile(ins.secondOnly, outfile, header);

    outfile = filenameBoth + "_Insertions.txt";
    header = "# Insertions in ";
    header += filenameBoth;
    utils::printVectorToFile(ins.both, outfile, header);

    //DELETIONS
    filenameEnding = "_OnlyDeletions.txt";
    outfile = filename + filenameEnding;
    header = "# Deletions in ";
    header += filename;
    header += " only";
    utils::printVectorToFile(dels.firstOnly, outfile, header);

    outfile = filename2 + filenameEnding;
    header = "# Deletions in ";
    header += filename2;
    header += " only";
    utils::printVectorToFile(dels.secondOnly, outfile, header);

    outfile = filenameBoth + "_Deletions.txt";
    header = "# Deletions in ";
    header += filenameBoth;
    utils::printVectorToFile(dels.both, outfile, header);

    //COMPLEX
    filenameEnding = "_OnlyComplex.txt";
    outfile = filename + filenameEnding;
    header = "# Complex in ";
    header += filename;
    header += " only";
    utils::printVectorToFile(comps.firstOnly, outfile, header);

    outfile = filename2 + filenameEnding;
    header = "# Complex in ";
    header += filename2;
    header += " only";
    utils::printVectorToFile(comps.secondOnly, outfile, header);

    outfile = filenameBoth + "_Complex.txt";
    header = "# Complex in ";
    header += filenameBoth;
    utils::printVectorToFile(comps.both, outfile, header);
  }
}

void CBasicStatistic::computeInsComparison(vector<InsRep> &aFirst, vector<InsRep> aSecond, CompareRes &aResult) {
  vector<InsRep> &firstOnly = mCompareINSRes.firstOnly;
  vector<InsRep> &secondOnly = mCompareINSRes.secondOnly;
  vector<InsRep> &bothFromFirst = mCompareINSRes.bothFromFirst;

  set_intersection(aFirst.begin(), aFirst.end(), aSecond.begin(), aSecond.end(), back_inserter(bothFromFirst), [this](const InsRep &aFirst, const InsRep &aSecond){ return this->compareInsRep(aFirst, aSecond);});

  set_intersection(aSecond.begin(), aSecond.end(), aFirst.begin(), aFirst.end(), back_inserter(mCompareINSRes.bothFromSecond), [this](const InsRep &aFirst, const InsRep &aSecond){ return this->compareInsRep(aFirst, aSecond);});

  set_difference(aFirst.begin(), aFirst.end(), bothFromFirst.begin(), bothFromFirst.end(), back_inserter(firstOnly), [this](const InsRep &aFirst, const InsRep &aSecond){ return this->compareInsRep(aFirst, aSecond);});

  set_difference(aSecond.begin(), aSecond.end(), bothFromFirst.begin(), bothFromFirst.end(), back_inserter(secondOnly), [this](const InsRep &aFirst, const InsRep &aSecond){ return this->compareInsRep(aFirst, aSecond);});

  for(auto itr = firstOnly.begin(); itr != firstOnly.end(); ++itr) {
    aResult.firstOnly.push_back((*itr).str());
  }
  for(auto itr = secondOnly.begin(); itr != secondOnly.end(); ++itr) {
    aResult.secondOnly.push_back((*itr).str());
  }
  for(auto itr = mCompareINSRes.bothFromFirst.begin(); itr != mCompareINSRes.bothFromFirst.end(); ++itr) {
    aResult.both.push_back((*itr).str());
  }
}

void CBasicStatistic::computeDelComparison(vector<DelRep> &aFirst, vector<DelRep> aSecond, CompareRes &aResult) {
  vector<DelRep> &firstOnly = mCompareDELRes.firstOnly;
  vector<DelRep> &secondOnly = mCompareDELRes.secondOnly;
  vector<DelRep> &bothFromFirst = mCompareDELRes.bothFromFirst;

  set_intersection(aFirst.begin(), aFirst.end(), aSecond.begin(), aSecond.end(), back_inserter(bothFromFirst), [this](const DelRep &aFirst, const DelRep &aSecond){ return this->compareDelRep(aFirst, aSecond);});

  set_intersection(aSecond.begin(), aSecond.end(), aFirst.begin(), aFirst.end(), back_inserter(mCompareDELRes.bothFromSecond), [this](const DelRep &aFirst, const DelRep &aSecond){ return this->compareDelRep(aFirst, aSecond);});

  set_difference(aFirst.begin(), aFirst.end(), bothFromFirst.begin(), bothFromFirst.end(), back_inserter(firstOnly), [this](const DelRep &aFirst, const DelRep &aSecond){ return this->compareDelRep(aFirst, aSecond);});

  set_difference(aSecond.begin(), aSecond.end(), bothFromFirst.begin(), bothFromFirst.end(), back_inserter(secondOnly), [this](const DelRep &aFirst, const DelRep &aSecond){ return this->compareDelRep(aFirst, aSecond);});

  for(auto itr = firstOnly.begin(); itr != firstOnly.end(); ++itr) {
    aResult.firstOnly.push_back((*itr).str());
  }
  for(auto itr = secondOnly.begin(); itr != secondOnly.end(); ++itr) {
    aResult.secondOnly.push_back((*itr).str());
  }
  for(auto itr = mCompareDELRes.bothFromFirst.begin(); itr != mCompareDELRes.bothFromFirst.end(); ++itr) {
    aResult.both.push_back((*itr).str());
  }
}

void CBasicStatistic::computeComparison(vector<string> &aFirst, vector<string> &aSecond, CompareRes &aResult) {
  set_intersection(aFirst.begin(), aFirst.end(), aSecond.begin(), aSecond.end(), back_inserter(aResult.both), [this](const string &aFirst, const string &aSecond){ return this->compareStringsWithoutGT(aFirst, aSecond);});

  set_difference(aFirst.begin(), aFirst.end(), aResult.both.begin(), aResult.both.end(), back_inserter(aResult.firstOnly), [this](const string &aFirst, const string &aSecond){ return this->compareStringsWithoutGT(aFirst, aSecond);});

  set_difference(aSecond.begin(), aSecond.end(), aResult.both.begin(), aResult.both.end(), back_inserter(aResult.secondOnly), [this](const string &aFirst, const string &aSecond){ return this->compareStringsWithoutGT(aFirst, aSecond);});
}

bool CBasicStatistic::compareStringsWithoutGT(const string &aFirst, const string &aSecond) {

  vector<string> firstEntries, secondEntries;
  boost::split(firstEntries, aFirst, boost::is_any_of("\t"));
  boost::split(secondEntries, aSecond, boost::is_any_of("\t"));

  int chrFirst = mOrderMap[firstEntries[0]];
  int chrSecond = mOrderMap[secondEntries[0]];

  if(chrFirst < chrSecond) {
    return true;
  } else if (chrFirst > chrSecond) {
    return false;
  }

  long long int posFirst = atoll(firstEntries[1].c_str());
  long long int posSecond = atoll(secondEntries[1].c_str());
  if(posFirst < posSecond) {
    return true;
  } else if (posFirst > posSecond) {
    return false;
  }
  return (firstEntries[2]).compare(secondEntries[2]);
}

bool CBasicStatistic::compareInsRep(const InsRep &aFirst, const InsRep &aSecond) {
  if(mOrderMap[aFirst.chrom] < mOrderMap[aSecond.chrom]) {
    return true;
  }

  unsigned int posDiff = abs(aFirst.pos - aSecond.pos);
  unsigned int lengthDiff = abs(aFirst.length - aSecond.length);
  if(mOrderMap[aFirst.chrom] == mOrderMap[aSecond.chrom]) {
    if(posDiff == 0 && lengthDiff == 0) {
      return false;
    }
    if(posDiff > mThresholds.maxInsPosDiff) {
      return aFirst.pos < aSecond.pos;
    }
    if(lengthDiff > mThresholds.maxInsLengthDiff) {
      if(posDiff == 0) {
        return aFirst.length < aSecond.length;
      } else {
        return aFirst.pos < aSecond.pos;
      }
    }
  }
  return false;
}

bool CBasicStatistic::compareDelRep(const DelRep &aFirst, const DelRep &aSecond) {
  if(mOrderMap[aFirst.chrom] < mOrderMap[aSecond.chrom]) {
    return true;
  }
  unsigned int posDiff = abs(aFirst.pos - aSecond.pos);
  unsigned int centerDiff = abs(aFirst.center - aSecond.center);
  unsigned int lengthDiff = abs(aFirst.length - aSecond.length);
  if(mOrderMap[aFirst.chrom] == mOrderMap[aSecond.chrom]) {
    if(centerDiff == 0 && lengthDiff == 0) {
      return false;
    }
    if(centerDiff > mThresholds.maxDelCenterDiff) {
      if(posDiff == 0) {
        return aFirst.length < aSecond.length;
      } else {
        return aFirst.pos < aSecond.pos;
      }
    }
    if(lengthDiff > mThresholds.maxInsLengthDiff) {
      if(posDiff == 0) {
        return aFirst.length < aSecond.length;
      } else {
        return aFirst.pos < aSecond.pos;
      }
    }
  }
  return false;
}

void CBasicStatistic::doBasicStats(vcflib::VariantCallFile &aVariantFile, const string &aFilename, CFilter *apFilter, bool aPrintOut) {

  BasicStats stats;

  Variant record(aVariantFile);
  stats.samples = aVariantFile.sampleNames.size();
  EVariantType alleleType = DEFAULT;
  while (aVariantFile.getNextVariant(record)) {//For all records in file
    ++stats.records;
    //Check info, effect, sample and predefined filters
    if(!apFilter->passesInfoFilter(record)) {
      continue;
    }

    if(!apFilter->passesEffectFilter(record)) {
      continue;
    }

    auto firstSample = record.sampleNames.begin();
    if(!apFilter->passesSampleFilter(record, *firstSample)) {
      continue;
    }

    alleleType = DEFAULT;
    for (vector<string>::iterator alts = record.alt.begin(); alts != record.alt.end(); ++alts) {
      string& alternate = *alts;
      EVariantType varType =  getVariantType(record.ref, alternate);
      if(alleleType == DEFAULT) {
        alleleType = varType;
      }
      if(alleleType != varType) {
        alleleType = COM;
        break;
      }
    }
    auto itr = record.samples.begin();
    auto sample = itr->second;
    string& genotype = sample["GT"].front();
    if(!apFilter->passesPredefinedFilters(alleleType, genotype)) {
      continue;
    }

    if(aPrintOut) {
      cout << record << endl;
    }

    ++stats.passedFilter;

    normalizeGenotype(genotype);
    string recordString = record.vrepr();
    stats.VrepToGT[recordString] = genotype;

    if(mOrderMap.count(record.sequenceName) == 0) {
      mOrderMap[record.sequenceName] = mOrderMap.size();
    }

    //Count genotype information
    itr = record.samples.begin();
    auto itrEnd = record.samples.end();
    for(; itr != itrEnd; ++itr) {
      auto sample = itr->second;
      string& genotype = sample["GT"].front();
      stats.genotypes[genotype]++;
    }

    //Do variant specific statistic
    if(record.alt.size() > 1) {
      ++stats.multiAll;
    }
    int allIndex = 1;
    switch(alleleType) {
      case SNP:
        stats.snpsVec.push_back(recordString);
        allIndex = 1;
        for (vector<string>::iterator alts = record.alt.begin(); alts != record.alt.end(); ++alts) {
          ESubstitutionType substType = getSubstitutionType(record.ref, *alts);
          double alleleFreq = computeAlleleFreq(record.samples, allIndex++);
          stats.substTypes[substType]++;
          double substFreq = stats.substFrequency[substType];
          if(substFreq < 0.0) {
            stats.substFrequency[substType] = alleleFreq;
          } else {
            stats.substFrequency[substType] = (stats.substFrequency[substType] + alleleFreq)/2.0;
          }
        }
        ++stats.snps;
        if(record.alt.size() > 1) {
          ++stats.multiAllSNPs;
        }
        break;
      case MNP:
        stats.mnpsVec.push_back(recordString);
        ++stats.mnps;
        break;
      case INS:
      {
        auto insertionRep = computeInsRep(record, genotype);
        stats.insRepVec.push_back(insertionRep);
        if(insertionRep.length > 50) {
          ++stats.largerIns;
        } else {
          updateIndelLengthFrequency(record, stats);
          ++stats.ins;
        }
        break;
      }
      case DEL:
      {
        auto deletionRep = computeDelRep(record, genotype);
        stats.delRepVec.push_back(deletionRep);

        if(deletionRep.length > 50) {
          ++stats.largerDels;
        } else {
          updateIndelLengthFrequency(record, stats);
          ++stats.dels;
        }
        break;
      }
      default:
        stats.compVec.push_back(recordString);
        ++stats.complexVar;
        break;
    }
  }
  mStatsPerFile[aFilename] = stats;
}

void CBasicStatistic::printGenotype(const string &aFilename) {
  BasicStats stats = mStatsPerFile[aFilename];
  auto itr = stats.genotypes.begin();
  auto itrEnd = stats.genotypes.end();
  for(; itr != itrEnd; ++itr) {
    cout << itr->first << " " << itr->second << endl;
  }
}

void CBasicStatistic::printBasicStats(const string &aFilename) {
  BasicStats stats = mStatsPerFile[aFilename];
  const char separator = ' ';
  const int nameWidth = 25;
  const int numWidth = 5;
  cout << left << setw(nameWidth) << setfill(separator) << "samples" << setw(numWidth) << setfill(separator) << stats.samples << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "total records" << setw(numWidth) << setfill(separator) << stats.records << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "records passed filter: " << setw(numWidth) << setfill(separator) << stats.passedFilter << " (" << fixed << setprecision(2) << ((double) stats.passedFilter/(double) stats.records)*100.0 << " % of total)" << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "snps" << setw(numWidth) << setfill(separator) << stats.snps << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "mnps" << setw(numWidth) << setfill(separator) << stats.mnps << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "insertions" << setw(numWidth) << setfill(separator) << (stats.ins + stats.largerIns) << " (including " << stats.largerIns << " larger insertions)" << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "deletions" << setw(numWidth) << setfill(separator) << (stats.dels + stats.largerDels) << " (including " << stats.largerDels << " larger deletions)" << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "complex" << setw(numWidth) << setfill(separator) << stats.complexVar << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "multiallelic sites" << setw(numWidth) << setfill(separator) << stats.multiAll << endl;
  cout << left << setw(nameWidth) << setfill(separator) << "multiallelic snp sites" << setw(numWidth) << setfill(separator) << stats.multiAllSNPs << endl;
}

map<string, string> CBasicStatistic::getVreprToGT(const string &aFilename) {
  return mStatsPerFile[aFilename].VrepToGT;
}

vector<pair<string, unsigned int>> CBasicStatistic::getSubstData(const string &aFilename) {
  BasicStats stats = mStatsPerFile[aFilename];
  bool substData = false;
  for(unsigned int i = AT; i < CT+1; ++i) {
    if(stats.substTypes[i] > 0) {
      substData = true;
    }
    break;
  }

  vector<pair<string, unsigned int>> aData;
  if(substData) {
    for(unsigned int i = AT; i < CT+1; ++i) {
     aData.push_back(make_pair(substNames[i], stats.substTypes[i]));
    }
  }
  return aData;
}

vector<pair<int, unsigned int>> CBasicStatistic::getIndelData(const string &aFilename) {
  BasicStats stats = mStatsPerFile[aFilename];
  vector<pair<int, unsigned int>> aData(stats.indelLengthCounts.size());
  for(auto itr = stats.indelLengthCounts.begin(); itr != stats.indelLengthCounts.end(); ++itr) {
    aData.push_back(make_pair(itr->first, itr->second));
  }
  return aData;
}

void CBasicStatistic::writeBasicStatsToFile(const string &aFilename, const string &aOutfile) {
  BasicStats stats = mStatsPerFile[aFilename];
  unsigned int snps = stats.snps;
  unsigned int mnps = stats.mnps;
  unsigned int ins = stats.ins;
  unsigned int dels = stats.dels;
  unsigned int complexVar = stats.complexVar;
  unsigned int data = snps + mnps + ins + dels + complexVar;
  if(data > 0) {
    std::ofstream outfile(aOutfile);
    if(!outfile.is_open()) {
      std::cerr << "Warning: could not open '" << aOutfile << "'" << std::endl;
      std::cerr << "Maybe some plots cannot be created or gnuplot will fail" << endl;
      return;
    }
    string filename = aFilename.substr(aFilename.find_last_of("\\/")+1);
    outfile << "file snp mnp insertion deletion complex" << endl;
    outfile << filename << " " << snps << " " << mnps << " " << ins << " " << dels << " " << complexVar;
    outfile.close();
  }
}

void CBasicStatistic::writeIndelFreqDataToFile(const string &aFilename, const string &aOutfile) {
  BasicStats stats = mStatsPerFile[aFilename];
  if(!stats.indelLengthCounts.empty() && !stats. indelLengthFrequency.empty()) {
    std::ofstream outfile(aOutfile);
    if(!outfile.is_open()) {
      std::cerr << "Warning: could not open '" << aOutfile << "'" << std::endl;
      std::cerr << "Maybe some plots cannot be created or gnuplot will fail" << endl;
      return;
    }
    string filename = aFilename.substr(aFilename.find_last_of("\\/")+1);
    outfile << "length counts frequency" << endl;
    for(auto itr = stats.indelLengthCounts.begin(); itr != stats.indelLengthCounts.end(); ++itr) {
      outfile << itr->first << " " << itr->second << " " << stats.indelLengthFrequency[itr->first] << endl;
    }
    outfile.close();
  }
}


void CBasicStatistic::writeSubstFreqDataToFile(const string &aFilename, const string &aOutfile) {
  BasicStats stats = mStatsPerFile[aFilename];
  bool hasTypeCounts = false;
  bool hasFreqs = false;
  for(unsigned int i = AT; i < CT+1; ++i) {
    if(stats.substTypes[i] != 0) {
      hasTypeCounts = true;
    }
    if(stats.substFrequency[i] != -1.0) {
      hasFreqs = true;
    }
  }
 if(hasTypeCounts && hasFreqs) {
    std::ofstream outfile(aOutfile);
    if(!outfile.is_open()) {
      std::cerr << "Warning: could not open '" << aOutfile << "'" << std::endl;
      std::cerr << "Maybe some plots cannot be created or gnuplot will fail" << endl;
      return;
    }
    string filename = aFilename.substr(aFilename.find_last_of("\\/")+1);
    outfile << "subst counts frequency" << endl;
    for(unsigned int i = AT; i < CT+1; ++i) {
      outfile << substNames[i] << " " << stats.substTypes[i] << " " << stats.substFrequency[i] << endl;
    }
    outfile.close();
  }
}

double CBasicStatistic::computeAlleleFreq(map<string, map<string, vector<string> > > aSamples, int aAllIndex) {
  double alleles = 0.0;
  double freq = 0.0;
  auto samplesItr = aSamples.begin();
  auto samplesItrEnd = aSamples.end();
  for(; samplesItr != samplesItrEnd; ++samplesItr) {
    alleles += 2.0;
    auto sample = samplesItr->second;
    string& genotype = sample["GT"].front();
    vector<string> gt = split(genotype, "|/");
    if(atoi(gt[0].c_str()) == aAllIndex) {
      freq++;
    }
    if(atoi(gt[1].c_str()) == aAllIndex) {
      freq++;
    }
  }
  return freq/alleles;
}

void CBasicStatistic::updateIndelLengthFrequency(Variant &aRecord, BasicStats &aBasicStats) {
  int allIndex = 1;
  for (vector<string>::iterator alts = aRecord.alt.begin(); alts != aRecord.alt.end(); ++alts) {
    int length = (*alts).size() - aRecord.ref.size();
    aBasicStats.indelLengthCounts[length]++;
    double alleleFreq = computeAlleleFreq(aRecord.samples, allIndex++);
    if(length <= 50 && length >= -50) {
      if(aBasicStats.indelLengthFrequency.count(length) == 0) {
        aBasicStats.indelLengthFrequency[length] = alleleFreq;
      } else {
        aBasicStats.indelLengthFrequency[length] = (aBasicStats.indelLengthFrequency[length] + alleleFreq)/2.0;
      }
    }
  }
}

void CBasicStatistic::normalizeGenotype(string &aGenotype) {
  vector<string> gt = split(aGenotype, "/");
  if(gt.size() == 2) {
    if(atoi(gt[1].c_str()) < atoi(gt[0].c_str())) {
      aGenotype = gt[1];
      aGenotype += "/";
      aGenotype += gt[0];
    }
  }
}

CBasicStatistic::InsRep CBasicStatistic::computeInsRep(Variant aVariant, const string &aGenotype) {
  InsRep insRep;
  insRep.chrom = aVariant.sequenceName;
  insRep.length = (*aVariant.alt.begin()).size() - aVariant.ref.size();
  insRep.pos = aVariant.position;
  insRep.genotype = aGenotype;
  return insRep;
}

CBasicStatistic::DelRep CBasicStatistic::computeDelRep(Variant aVariant, const string &aGenotype) {
  DelRep delRep;
  delRep.chrom = aVariant.sequenceName;
  delRep.length = aVariant.ref.size() - (*aVariant.alt.begin()).size();
  delRep.center = aVariant.position + ((int) delRep.length / 2) + 1;
  delRep.pos = aVariant.position;
  delRep.genotype = aGenotype;
  return delRep;
}
