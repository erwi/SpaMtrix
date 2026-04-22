#include <catch.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "spamtrix_ircmatrix.hpp"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_reader.hpp"
#include "spamtrix_vector.hpp"
#include "spamtrix_writer.hpp"

namespace {
using namespace SpaMtrix;

std::filesystem::path makeTempPath(const std::string &stem, const std::string &extension) {
  return std::filesystem::temp_directory_path() / (stem + extension);
}

void removeFileIfExists(const std::filesystem::path &path) {
  std::error_code ec;
  std::filesystem::remove(path, ec);
}

void requireMatricesEqual(const IRCMatrix &expected, const IRCMatrix &actual) {
  REQUIRE(expected.getNumRows() == actual.getNumRows());
  REQUIRE(expected.getNumCols() == actual.getNumCols());
  REQUIRE(expected.getnnz() == actual.getnnz());

  for (idx row = 0; row < expected.getNumRows(); ++row) {
    for (idx col = 0; col < expected.getNumCols(); ++col) {
      REQUIRE(expected.getValue(row, col) == Approx(actual.getValue(row, col)));
    }
  }
}

Vector readCsvVector(const std::filesystem::path &path) {
  std::ifstream input(path);
  REQUIRE(input.is_open());

  std::vector<real> values;
  std::string line;
  while (std::getline(input, line)) {
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
      if (!token.empty()) {
        values.push_back(std::stod(token));
      }
    }
  }

  Vector result(values.size());
  for (idx i = 0; i < values.size(); ++i) {
    result[i] = values[i];
  }
  return result;
}
}

TEST_CASE("Reader parses a MatrixMarket fixture", "[Reader][io]") {
  using namespace SpaMtrix;

  const auto path = makeTempPath("spamtrix_reader_fixture", ".mtx");
  removeFileIfExists(path);

  {
    std::ofstream output(path);
    REQUIRE(output.is_open());
    output << "% MatrixMarket matrix coordinate real general\n";
    output << "% tiny fixture\n";
    output << "2 3 2\n";
    output << "1 2 4.5\n";
    output << "2 3 -1.25\n";
  }

  IRCMatrix matrix = Reader::readMatrixMarket(path.string());
  REQUIRE(matrix.getNumRows() == 2);
  REQUIRE(matrix.getNumCols() == 3);
  REQUIRE(matrix.getnnz() == 2);
  REQUIRE(matrix.getValue(0, 1) == Approx(4.5));
  REQUIRE(matrix.getValue(1, 2) == Approx(-1.25));
  REQUIRE(matrix.getValue(0, 0) == Approx(0.0));

  removeFileIfExists(path);
}

TEST_CASE("MatrixMarket writer and reader preserve matrix contents", "[Reader][Writer][io]") {
  using namespace SpaMtrix;

  MatrixMaker mm(3, 4);
  mm.addNonZero(0, 0, 1.0);
  mm.addNonZero(0, 3, 2.0);
  mm.addNonZero(1, 1, 3.0);
  mm.addNonZero(2, 2, 4.0);
  IRCMatrix expected = mm.getIRCMatrix();

  const auto path = makeTempPath("spamtrix_matrix_roundtrip", ".mtx");
  removeFileIfExists(path);

  REQUIRE(Writer::writeMatrixMarket(path.string(), expected));
  IRCMatrix actual = Reader::readMatrixMarket(path.string());
  requireMatricesEqual(expected, actual);

  removeFileIfExists(path);
}

TEST_CASE("CSV writer produces parseable output", "[Writer][io]") {
  using namespace SpaMtrix;

  Vector expected(6);
  for (idx i = 0; i < expected.getLength(); ++i) {
    expected[i] = static_cast<real>(i) * 1.5;
  }

  const auto path = makeTempPath("spamtrix_vector_roundtrip", ".csv");
  removeFileIfExists(path);

  REQUIRE(Writer::writeCSV(path.string(), expected, 3));
  Vector actual = readCsvVector(path);

  REQUIRE(actual.getLength() == expected.getLength());
  for (idx i = 0; i < expected.getLength(); ++i) {
    REQUIRE(actual[i] == Approx(expected[i]));
  }

  removeFileIfExists(path);
}