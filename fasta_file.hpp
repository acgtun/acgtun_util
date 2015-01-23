#ifndef FASTAFILE_H_
#define FASTAFILE_H_

#include "sdk.hpp"
#include "bio_util.hpp"

#include <vector>
using std::vector;

class FastaFile {
 public:
  FastaFile(const string& file_path);
  void SetFilePath(const string& file_path);
  string GetFilePath() const;
  FastaFile();
  ~FastaFile();

 private:
  uint64_t GetLineFromString(const char* strVal, char* strRet);
  void GetNumOfSequences();
  void AnalyzeFile();
  void ReadFastaFile();

  string file_path;  // Absolute path for the fasta file
  char* file_string;  // To store characters of the whole file
  uint64_t file_size;

 public:
  vector<string> sequences_names;
  char** sequences;
  uint32_t num_of_sequences;
  uint64_t num_of_characters;
  uint32_t max_sequence_length;
};

#endif /* FASTAFILE_H_ */
