#include <boost/test/unit_test.hpp>
#include "../src/Filter.hpp"

using namespace std;
using namespace vcflib;

BOOST_AUTO_TEST_CASE( FilterTest ) {
  map<string, VariantFieldType> infoMap, sampleMap;
  string infoFilter = "";
  string sampleFilter = "";
  string effectFilter = "";
  string sample = "";
  CFilter::PredefinedFilters predefFilters;
  CFilter *filterT = new CFilter(infoFilter, infoMap, sampleFilter, sampleMap, effectFilter, predefFilters);
  Variant var;
  BOOST_CHECK_EQUAL(filterT->passesInfoFilter(var), true);
  BOOST_CHECK_EQUAL(filterT->passesSampleFilter(var, sample), true);
  delete filterT;

  vcflib::VariantCallFile vcffile;
  string filename = "testdata/short.vcf";
  bool error = false;
  try{
    vcffile.open(filename);
  }
  catch(...) {
    error = true;
  }
  BOOST_CHECK_EQUAL(false, error);

  Variant record(vcffile);
  vcffile.getNextVariant(record);

  infoFilter = "Invalid = 2";
  sampleFilter = "Illegal < 9";

  streambuf* oldCerr = cerr.rdbuf();
  cerr.rdbuf(0);
  CFilter *filterF = new CFilter(infoFilter, vcffile.infoTypes, sampleFilter, vcffile.formatTypes, effectFilter, predefFilters);
  BOOST_CHECK_EQUAL(filterF->passesInfoFilter(record), false);
  BOOST_CHECK_EQUAL(filterF->passesSampleFilter(record, sample), false);
  cerr.rdbuf(oldCerr);
  delete filterF;

  infoFilter = "DP = 14";
  sampleFilter = "GT = 0|0";
  sample = "NA00001";
  CFilter *filter = new CFilter(infoFilter, vcffile.infoTypes, sampleFilter, vcffile.formatTypes, effectFilter, predefFilters);
  BOOST_CHECK_EQUAL(filter->passesInfoFilter(record), true);
  BOOST_CHECK_EQUAL(filter->passesSampleFilter(record, sample), true);
  sample = "NA00002";
  BOOST_CHECK_EQUAL(filter->passesSampleFilter(record, sample), false);
  delete filter;
}
