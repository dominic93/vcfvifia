#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "../src/utils.hpp"
#include <fstream>
#include <sstream>
#include <string>


using namespace std;
namespace bf = boost::filesystem;

BOOST_AUTO_TEST_CASE( UtilsTest_printVectorToFile ) {
  bf::path tempDir = bf::unique_path();
  bf::path tempFile(tempDir /"temp.dat");
  bf::path tempFile2(tempDir /"temp2.dat");
  BOOST_CHECK_EQUAL(bf::create_directory(tempDir), true);
  vector<int> intVec;
  BOOST_CHECK_EQUAL(utils::printVectorToFile(intVec, tempFile.c_str()), true);
  BOOST_CHECK_EQUAL(bf::exists(tempFile), false);

  intVec.push_back(1);
  BOOST_CHECK_EQUAL(utils::printVectorToFile(intVec, tempFile.c_str()), true);
  BOOST_CHECK_EQUAL(bf::exists(tempFile), true);

  vector<string> stringVec;
  stringVec.push_back("42");
  stringVec.push_back("Test strings");
  stringVec.push_back("?");
  BOOST_CHECK_EQUAL(utils::printVectorToFile(stringVec, tempFile2.c_str()), true);
  BOOST_CHECK_EQUAL(bf::exists(tempFile2), true);

  ifstream infile(tempFile2.c_str());
  string line;
  int lineNr = 0;
  while (getline(infile, line)) {
    if(lineNr == 3) {
      if(line != "") {
        BOOST_CHECK_EQUAL(false, lineNr);
      }
    } else {
      BOOST_CHECK_EQUAL(line, stringVec.at(lineNr++));
    }
  }
  BOOST_CHECK_EQUAL(lineNr, 3);

  infile.close();
  remove_all(tempDir);
}

BOOST_AUTO_TEST_CASE( UtilsTest_firstIsDigit ) {
  string firstIs = "1asdsad";
  string firstIsnt = "i";
  string empty = "";
  BOOST_CHECK_EQUAL(true, utils::firstIsDigit(firstIs));
  BOOST_CHECK_EQUAL(false, utils::firstIsDigit(firstIsnt));
  BOOST_CHECK_EQUAL(false, utils::firstIsDigit(empty));
}
