#ifndef BIO_UTIL_H_
#define BIO_UTIL_H_

#include <string>

using namespace std;

/*
 A0
 R1
 N2
 D3
 C4
 E5
 Q6
 G7
 H8
 I9
 L10
 K11
 M12
 F13
 P14
 S15
 T16
 W17
 Y18
 V19
 */
/* B J O U X Z*/

#define HASHLEN 6
const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";
                   /*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */
const int base[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };
const usint32_t basep[] = { 1, 20, 400, 8000, 160000, 3200000, 64000000, 1280000000 };

const int BLOSUM62[][20] = {
//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 },  //A
{ -1, 5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 },  //R
{ -2, 0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 },  //N
{ -2,-2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 },  //D
{  0,-3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },  //C
{ -1, 1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 },  //Q
{ -1, 0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 },  //E
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 },  //G
{ -2, 0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 },  //H
{ -1,-3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 },  //I
{ -1,-2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 },  //L
{ -1, 2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 },  //K
{ -1,-1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 },  //M
{ -2,-3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 },  //F
{ -1,-2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 },  //P
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 },  //S
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 },  //T
{ -3,-3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 },  //W
{ -2,-2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 },  //Y
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3,  -1, 4 } };  //V
//  http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#get_subsequence
//  http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

/* Given a 6-mer amino acids sequence, translate it to a integer number.
 * Use 20-based number, A is 0, R is 1, N is 3 and so on.
 * */
inline uint32_t GetHashValueAA(const char * seed) {
  usint32_t hash_value = 0;
  for (usint32_t i = 0; i < HASHLEN; i++) {
    hash_value += base_preCal[base[seed[i] - 'A']][i];
  }
  return hash_value;
}

/* Given a integer, translate it to a 6-mer amino acids, it's also 20-based.*/
void DeCodeAA(const usint32_t& hash_value, string& seed) {
  usint32_t value = hash_value;
  for (usint32_t i = 0; i < HASHAALEN; i++) {
    seed += AA20[value % 20];
    value /= 20;
  }
}

#endif /* BIO_UTIL_H_ */
