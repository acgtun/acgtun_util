#include "fasta_file.hpp"

void FastaFile::SetFilePath(const string& _file_path) {
  file_path = _file_path;
  ReadFastaFile();
}

string FastaFile::GetFilePath() const {
  return file_path;
}

void FastaFile::ReadFastaFile() {
  INFO("Read Fasta file:", file_path.c_str());
  FILE* fin = fopen(file_path.c_str(), "rb");
  FILE_OPEN_CHECK(fin);
  fseek(fin, 0, SEEK_END);
  file_size = ftell(fin);
  INFO("file_size", file_size);
  MEMORY_ALLOCATE_CHECK(
      file_string = (char*) malloc(sizeof(char) * (file_size + 1)));
  fseek(fin, 0, SEEK_SET);
  FREAD_CHECK(fread(file_string, 1, file_size, fin), file_size);
  fclose(fin);

  sequences = NULL;
  num_of_sequences = 0;
  num_of_characters = 0;

  GetNumOfSequences();
  AnalyzeFile();
}

FastaFile::FastaFile(const string& _file_path) {
  file_path = _file_path;
  file_string = NULL;
  file_size = 0;

  sequences = NULL;
  num_of_sequences = 0;
  num_of_characters = 0;
  max_sequence_length = 100000;
  ReadFastaFile();
}

FastaFile::FastaFile() {
  file_path = "";
  file_string = NULL;
  file_size = 0;

  sequences = NULL;
  num_of_sequences = 0;
  num_of_characters = 0;
  max_sequence_length = 100000;
}

FastaFile::~FastaFile() {
  for (uint32_t i = 0; i < num_of_sequences; i++) {
    free (sequences[i]);
  }
  free (sequences);
  sequences_names.clear();
}

uint64_t FastaFile::GetLineFromString(const char* strVal, char* strRet) {
  uint64_t i;
  bool tag = 0;
  uint64_t j = 0;
  for (i = 0; strVal[i] != 0; i++) {
    if (strVal[i] == ' ')
      tag = 1;
    if (0xA == strVal[i] || 0xD == strVal[i]) {
      break;
    }
    if (tag == 0) {
      strRet[j] = strVal[i];
      j++;
    }
  }

  strRet[j] = 0;
  return i;
}

void FastaFile::GetNumOfSequences() {
  char strRet[MAX_LINE_LEN];
  for (uint64_t i = 0; i < file_size; i++) {
    if (file_string[i] == '>') {
      num_of_sequences++;
    }
    i += GetLineFromString(&file_string[i], strRet);
  }
}

void FastaFile::AnalyzeFile() {
  MEMORY_ALLOCATE_CHECK(
      sequences = (char**) malloc(sizeof(char *) * num_of_sequences));
  sequences_names.resize(num_of_sequences);

  uint32_t id = 0;
  string str_seq;
  num_of_characters = 0;
  max_sequence_length = 0;
  str_seq.clear();
  char strRet[MAX_LINE_LEN];
  for (uint64_t i = 0; i < file_size; i++) {
    if (file_string[i] == '>') {
      if (str_seq.size() != 0) {
        MEMORY_ALLOCATE_CHECK(
            sequences[id] = (char *) malloc(
                sizeof(char) * (str_seq.size() + 1)));
        strcpy(sequences[id], str_seq.c_str());
        num_of_characters += str_seq.size();
        if (str_seq.size() > max_sequence_length) {
          max_sequence_length = str_seq.size();
        }
        str_seq.clear();
        id++;
      }
      i += GetLineFromString(&file_string[i], strRet);
      sequences_names[id] = string(strRet).substr(1);
    } else if (AA20.find_first_of(file_string[i]) != string::npos) {
      str_seq += toupper(file_string[i]);
    } else if (isalpha (file_string[i])) {
      str_seq += AA20[rand() % 20];  //todo
    }
  }
  if (str_seq.size() != 0) {
    MEMORY_ALLOCATE_CHECK(
        sequences[id] = (char *) malloc(sizeof(char) * (str_seq.size() + 1)));
    strcpy(sequences[id], str_seq.c_str());
    num_of_characters += str_seq.size();
    if (str_seq.size() > max_sequence_length) {
      max_sequence_length = str_seq.size();
    }
    str_seq.clear();
    id++;
  }
  if (id != num_of_sequences) {
    ERROR_INFO("The number of Sequences has error~!");
  }

  free (file_string);

  //////////////////////////////////////////////////////
  FILE * fout = fopen("proteins.txt", "w");
  for (uint32_t i = 0; i < num_of_sequences; i++) {
    fprintf(fout, "%d:%s\n", i, sequences_names[i].c_str());
    fprintf(fout, "%s\n", sequences[i]);
  }
  fclose(fout);
  ///////////////////////////////////////////////////////
}
