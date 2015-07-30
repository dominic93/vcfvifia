#include <boost/test/unit_test.hpp>
#include "../src/BasicStatistic.hpp"
#include "utils.hpp"
#include <boost/filesystem.hpp>

using namespace std;
namespace bf = boost::filesystem;

BOOST_AUTO_TEST_CASE( BasicStatisticTest_doAndPrintBasicStats ) {
  CBasicStatistic::Thresholds thresholds;
  CBasicStatistic basicStatistic(thresholds);
  string filename = "testdata/short.vcf";
  vcflib::VariantCallFile file;
  file.open(filename);
  CFilter::PredefinedFilters predefFilters;
  CFilter *filter;
  filter = new CFilter(predefFilters);
  basicStatistic.doBasicStats(file, filename, filter, false);
  delete filter;

  bf::path tempDir = bf::unique_path();
  bf::path tempFile(tempDir/"temp.txt");
  BOOST_CHECK_EQUAL(bf::create_directory(tempDir), true);

  streambuf* oldCoutStreamBuf = cout.rdbuf();
  ofstream result(tempFile.c_str());
  cout.rdbuf( result.rdbuf() );
  basicStatistic.printBasicStats(filename);
  cout.rdbuf( oldCoutStreamBuf );

  string filenameExp = "testdata/expResults/output1.txt";
  ifstream expected(filenameExp);
  ifstream comp(tempFile.c_str());
  istream_iterator<char> begin(expected);
  istream_iterator<char> begin2(comp);
  istream_iterator<char> end;
  BOOST_CHECK_EQUAL_COLLECTIONS(begin, end, begin2, end);

  result.close();
  expected.close();
  comp.close();
  remove_all(tempDir);
}

BOOST_AUTO_TEST_CASE( BasicStatisticTest_compareFiles ) {
  CBasicStatistic::Thresholds thresholds;
  CBasicStatistic basicStatistic(thresholds);
  string filename = "testdata/short.vcf";
  vcflib::VariantCallFile file;
  file.open(filename);
  CFilter::PredefinedFilters predefFilters;
  CFilter *filter;
  filter = new CFilter(predefFilters);
  basicStatistic.doBasicStats(file, filename, filter, false);
  delete filter;

  streambuf* oldCoutStreamBuf = cout.rdbuf();
  cout.rdbuf(0);
  basicStatistic.compareTwo(filename, filename, false, false);
  cout.rdbuf( oldCoutStreamBuf );
  CBasicStatistic::CompareRes compRes = basicStatistic.getCompareResByType(SNP);
  BOOST_CHECK_EQUAL(0, compRes.firstOnly.size());
  BOOST_CHECK_EQUAL(0, compRes.secondOnly.size());
  BOOST_CHECK_EQUAL(3, compRes.both.size());
}

BOOST_AUTO_TEST_CASE( BasicStatisticTest_getter ) {
  CBasicStatistic::Thresholds thresholds;
  CBasicStatistic basicStatistic(thresholds);
  string filename = "testdata/short.vcf";
  vcflib::VariantCallFile file;
  file.open(filename);
  CFilter::PredefinedFilters predefFilters;
  CFilter *filter;
  filter = new CFilter(predefFilters);
  basicStatistic.doBasicStats(file, filename, filter, false);
  delete filter;

  auto substData = basicStatistic.getSubstData(filename);
  auto indelData = basicStatistic.getIndelData(filename);
  int substCount = 0;
  int indelCount = 0;
  for(auto itr = substData.begin(); itr != substData.end(); ++itr) {
    substCount += (*itr).second;
  }
  for(auto itr = indelData.begin(); itr != indelData.end(); ++itr) {
    indelCount += (*itr).second;
  }
  BOOST_CHECK_EQUAL(4, substCount);
  BOOST_CHECK_EQUAL(2, indelCount);
}

BOOST_AUTO_TEST_CASE( BasicStatisticTest_writeStatisticToFile ) {
  CBasicStatistic::Thresholds thresholds;
  CBasicStatistic basicStatistic(thresholds);
  string filename = "testdata/short.vcf";
  vcflib::VariantCallFile file;
  file.open(filename);
  CFilter::PredefinedFilters predefFilters;
  CFilter *filter;
  filter = new CFilter(predefFilters);
  basicStatistic.doBasicStats(file, filename, filter, false);
  delete filter;

  bf::path tempDir = bf::unique_path();
  bf::path tempFile(tempDir /"temp.txt");
  bf::path tempFile2(tempDir /"temp1.txt");
  bf::path tempFile3(tempDir /"temp2.txt");

  BOOST_CHECK_EQUAL(bf::create_directory(tempDir), true);


  string filenameExp = "testdata/expResults/output2.txt";
  string filenameExp1 = "testdata/expResults/output3.txt";
  string filenameExp2 = "testdata/expResults/output4.txt";


  basicStatistic.writeBasicStatsToFile(filename, tempFile.c_str());
  basicStatistic.writeSubstFreqDataToFile(filename, tempFile2.c_str());
  basicStatistic.writeIndelFreqDataToFile(filename, tempFile3.c_str());

  ifstream expected(filenameExp);
  ifstream expected1(filenameExp1);
  ifstream expected2(filenameExp2);

  ifstream comp(tempFile.c_str());
  ifstream comp1(tempFile2.c_str());
  ifstream comp2(tempFile3.c_str());

  istream_iterator<char> begin0(expected);
  istream_iterator<char> begin02(comp);
  istream_iterator<char> end0;

  istream_iterator<char> begin1(expected1);
  istream_iterator<char> begin12(comp1);
  istream_iterator<char> end1;

  istream_iterator<char> begin2(expected2);
  istream_iterator<char> begin22(comp2);
  istream_iterator<char> end2;


  BOOST_CHECK_EQUAL_COLLECTIONS(begin0, end0, begin02, end0);
  BOOST_CHECK_EQUAL_COLLECTIONS(begin1, end1, begin12, end1);
  BOOST_CHECK_EQUAL_COLLECTIONS(begin2, end2, begin22, end2);

  expected.close();
  expected1.close();
  expected2.close();
  comp.close();
  comp1.close();
  comp2.close();
  remove_all(tempDir);
}
