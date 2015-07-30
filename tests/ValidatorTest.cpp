#include <boost/test/unit_test.hpp>
#include "../src/BasicStatistic.hpp"
#include "../src/Validator.hpp"
#include "utils.hpp"
#include <boost/filesystem.hpp>

using namespace std;
namespace bf = boost::filesystem;

BOOST_AUTO_TEST_CASE( ValidatorTest_computeAndPrintValdiationResults ) {
  CBasicStatistic::Thresholds thresholds;
  CBasicStatistic basicStatistic(thresholds);
  string filename = "testdata/venter.chr1.unphased.vcf";
  string filename2 = "testdata/freebayes.len50000.vcf";
  vcflib::VariantCallFile file, file2;
  file.open(filename);
  file2.open(filename2);
  CFilter::PredefinedFilters predefFilters;
  CFilter *filter;
  filter = new CFilter(predefFilters);
  basicStatistic.doBasicStats(file, filename, filter, false);
  basicStatistic.doBasicStats(file2, filename2, filter, false);

  bf::path tempDir = bf::unique_path();
  bf::path tempFile(tempDir/"temp.txt");
  BOOST_CHECK_EQUAL(bf::create_directory(tempDir), true);

  streambuf *oldCoutStreamBuf = cout.rdbuf();
  ofstream result(tempFile.c_str());
  cout.rdbuf( result.rdbuf() );

  basicStatistic.compareTwo(filename, filename2, false, false);
  delete filter;

  CValidator validator;
  vector<EVariantType> variants = {SNP, MNP, COM};
  for(auto itr = variants.begin(); itr != variants.end(); ++itr) {
    validator.computeValidationResults(basicStatistic.getCompareResByType(*itr), *itr, basicStatistic.getVreprToGT(filename), basicStatistic.getVreprToGT(filename2));
  }
  validator.computeValidationResults(basicStatistic.getCompareINSRes());
  validator.computeValidationResults(basicStatistic.getCompareDELRes());

  validator.printValidationResults(filename, filename2);
  cout.rdbuf( oldCoutStreamBuf );

  string filenameExp = "testdata/expResults/output5.txt";
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

