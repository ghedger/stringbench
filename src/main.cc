/* MIT License
 *
 * Copyright (c) 2018 Greg Hedger
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <algorithm>
#include "stringbench_common.h"
#include "stringbench_flags.h"
#include "stringbench_log.h"

namespace stringbench {

// OutputPreamble
void OutputPreamble()
{
  using namespace std;
  cout << "StringBench" << endl;
  cout << "Copyright (c) 2018 Greg Hedger" << endl;
  cout << "MIT License" << endl;
  cout << endl;
  cout << "Permission is hereby granted, free of charge, to any person obtaining a copy" << endl;
  cout << "of this software and associated documentation files (the \"Software\"), to deal" << endl;
  cout << endl;
  cout << "in the Software without restriction, including without limitation the rights" << endl;
  cout << "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell" << endl;
  cout << "copies of the Software, and to permit persons to whom the Software is" << endl;
  cout << "furnished to do so, subject to the following conditions:" << endl;
  cout << endl;
  cout << "The above copyright notice and this permission notice shall be included in all" << endl;
  cout << "copies or substantial portions of the Software." << endl;
  cout << "Copyright (C) 2019 Gregory P. Hedger" << endl;
  cout << endl;
}

// PrintUsage
// Print basic program usage/invocation
void PrintUsage()
{
  using namespace std;
  OutputPreamble();
  cout << "" << endl;
  cout << "Usage:" << endl;
  cout << "\tstringbench [flags] word_a word_b" << endl;
}


//
// Distance helpers
//

inline int Min3(int a, int b, int c)
{
  if (a < b)
    if (a < c)
      return a;
    else
      return c;
  else
    if (b < c)
      return b;
    else
      return c;
}

inline int Min4(int a, int b, int c, int d)
{
  if (a < b)
    return Min3(a,c,d);
  else
    return Min3(b,c,d);
}

inline int Min(int a, int b)
{
  if (a < b)
    return a;
  return b;
}

inline int Max(int a, int b)
{
  if (a > b)
    return a;
  return b;
}


//
// Distance implementations
//

// CalcLevenshtein
// Calculates Levenshtein string edit distance -
//   insertions + deletions + substitutions
// It returns how "different" two strings are, effectively performing
// a commutative subtraction operation.
// Entry: string a
//        string b
// Exit:  Levenshtein distance
int CalcLevenshtein(const char *s1, const char *s2)
{
  int score = -1;
  if (s1 && s2) {
    int len1, len2, x, y, lastdiag, olddiag;
    len1 = strlen(s1);
    len2 = strlen(s2);
    int column[ len1 + 1 ];
    for (y = 1; y <= len1; y++)
      column[ y ] = y;
    for (x = 1; x <= len2; x++) {
      column[ 0 ] = x;
      for (y = 1, lastdiag = x - 1; y <= len1; y++) {
        olddiag = column[ y ];
        column[y] = Min3(column[ y ] + 1,
          column[ y - 1 ] + 1,
          lastdiag +
          (
           s1[ y - 1 ] == s2[ x - 1 ] ? 0 : 1
          )
        );
        lastdiag = olddiag;
      }
    }
    score = column[ len1 ];
  }
  return score;
}

// CalcDamerauLevenshtein
// Calculates the Damerau-Levenshtein string edit distance -
//   insertions + deletions + substitutions + transpositions
// Entry: string a
//        string b
// Exit: Damerau-Levenshtein distance
int CalcDamerauLevenshtein(const char *s1, const char* s2) {
#define DAMERAU_LEVENSHTEIN_ACCESS(i,j) dd[(i) * (m + 2) + (j)]
    int *dd;
    int i, j, cost, i1, j1, db;
    int m = strlen(s1);
    int n = strlen(s2);
    int infinity = n + m;
    int da[256 * sizeof(int)];    // size of ascii charset domain
    memset(da, 0, sizeof(da));

    if (!(dd = (int *) calloc((n + 2) * (m + 2), sizeof(int)))) {
      return -1;
    }

    DAMERAU_LEVENSHTEIN_ACCESS(0, 0) = infinity;
    for(i = 0; i < n + 1; ++i) {
      DAMERAU_LEVENSHTEIN_ACCESS(i + 1, 1) = i;
      DAMERAU_LEVENSHTEIN_ACCESS(i + 1, 0) = infinity;
    }

    for(j = 0; j < m + 1; ++j) {
      DAMERAU_LEVENSHTEIN_ACCESS(1, j + 1) = j;
      DAMERAU_LEVENSHTEIN_ACCESS(0, j + 1) = infinity;
    }

    for(i = 1; i < n + 1; ++i) {
      db = 0;
      for(j = 1; j < m + 1; ++j) {
        i1 = da[s2[j - 1]];
        j1 = db;
        cost = ((s1[i - 1] == s2[j - 1]) ? 0 : 1);
        if(!cost)
          db = j;
        {
          int a = dd[i*(m+2)+j] + cost;
          int b = DAMERAU_LEVENSHTEIN_ACCESS(i + 1, j) + 1;
          int c = DAMERAU_LEVENSHTEIN_ACCESS(i, j + 1) + 1;
          int d = DAMERAU_LEVENSHTEIN_ACCESS(i1, j1) + (i - i1 - 1) + 1 + (j - j1 - 1);
        DAMERAU_LEVENSHTEIN_ACCESS(i + 1, j + 1) =
          Min4(DAMERAU_LEVENSHTEIN_ACCESS(i, j) + cost,
              DAMERAU_LEVENSHTEIN_ACCESS(i + 1, j) + 1,
              DAMERAU_LEVENSHTEIN_ACCESS(i, j + 1) + 1,
              DAMERAU_LEVENSHTEIN_ACCESS(i1, j1) + (i - i1 - 1) + 1 + (j - j1 - 1));
        }
      }
      if ( m > i - 1) {
        da[s1[i - 1]] = i;
      }
    }
    cost = DAMERAU_LEVENSHTEIN_ACCESS(n + 1, m + 1);
    free(dd);
    return cost;
#undef DAMERAU_LEVENSHTEIN_ACCESS
}

// JaroSimilarity
// Calculates Jaro similarity
// Entry: string a
//        string b
// Exit: Jaro similarity
double CalcJaroSimilarity(const char *s1, const char *s2) {
    // length of the strings
    int str1_len = strlen(s1);
    int str2_len = strlen(s2);

    // if both strings are empty return 1
    // if only one of the strings is empty return 0
    if (str1_len == 0) return str2_len == 0 ? 1.0 : 0.0;

    // max distance between two chars to be considered matching
    // floor() is ommitted due to integer division rules
    int match_distance = (int) Max(str1_len, str2_len)/2 - 1;

    // arrays of bools that signify if that char in the matching string has a match
    bool *str1_matches = (bool *) calloc(str1_len, sizeof(bool));
    bool *str2_matches = (bool *) calloc(str2_len, sizeof(bool));

    // number of matches and transpositions
    double matches = 0.0;
    double transpositions = 0.0;

    // find the matches
    for (int i = 0; i < str1_len; i++) {
        // start and end take into account the match distance
        int start = Max(0, i - match_distance);
        int end = Min(i + match_distance + 1, str2_len);

        for (int k = start; k < end; k++) {
            // if s2 already has a match continue
            if (str2_matches[k]) continue;
            // if s1 and s2 are not
            if (s1[i] != s2[k]) continue;
            // otherwise assume there is a match
            str1_matches[i] = true;
            str2_matches[k] = true;
            matches++;
            break;
        }
    }

    // if there are no matches return 0
    if (matches == 0) {
        free(str1_matches);
        free(str2_matches);
        return 0.0;
    }

    // count transpositions
    int k = 0;
    for (int i = 0; i < str1_len; i++) {
        // if there are no matches in s1 continue
        if (!str1_matches[i]) continue;
        // while there is no match in s2 increment k
        while (!str2_matches[k]) k++;
        // increment transpositions
        if (s1[i] != s2[k]) transpositions++;
        k++;
    }

    // divide the number of transpositions by two as per the algorithm specs
    // this division is valid because the counted transpositions include both
    // instances of the transposed characters.
    transpositions /= 2.0;

    // free the allocated memory
    free(str1_matches);
    free(str2_matches);

    // return the Jaro distance
    return ((matches / str1_len) +
        (matches / str2_len) +
        ((matches - transpositions) / matches)) / 3.0;
}

// CalcJaroWinklerSimilarity
// The Jaro-Winkler similarity extense the Jaro similarity by giving
// additional weight to a match of the first (up to) four characters.
// Entry: string a
//        string b
// Exit: Jaro-Winkler similarity
double CalcJaroWinklerSimilarity(const char *s1, const char *s2)
{
  double scale = 0.1; // Winkler's original work used scaling factor 0.1
  double jaro_similarity = CalcJaroSimilarity(s1, s2);
  double jaro_winkler;

  int i = 0;
  while( i < 4 && s1[i] == s2[i]) {
    i++;
  }
  jaro_winkler = jaro_similarity + ((double) i * scale * (1 - jaro_similarity));
  return jaro_winkler;
}


// Cleanstring
// RTrim and lowercase a string
// Entry: string
void CleanString(std::string& word)
{
  using namespace std;
  word.erase(find_if(word.rbegin(), word.rend(), [](int c) {
        return !isspace(c);
    }).base(), word.end());
  transform(word.begin(), word.end(), word.begin(), ::tolower);
}
} // namespace stringbench

// main
// This is the main entry point and testbed for ternary tree.
//
// @In:     -
// @Out:    0 == success
int main(int argc, const char *argv[])
{
  using namespace std;
  using namespace stringbench;

  // Sets a sane level that allows UI but not debug messages
  SET_VERBOSITY_LEVEL(LOG_NORMAL);

  // This parses the arguments and takes subsequent non-dashed arguments
  // as the input (no quotes required)
  StringBenchFlags flags;
  flags.word = 0;
  //flags.bit.tree_engine = 0;
  //flags.bit.allow_dupes = flags.bit.output_directly
  //  = flags.bit.big_dictionary = 0;
  string word_a, word_b;
  if (1 < argc) {
    int i = 1;
    while (i < argc) {
      if ('-' == argv[i][0]) {
        switch(argv[i][1]) {
          case 'v': {
              int verbosity;
              if (isdigit(verbosity = argv[i][2])) {
                LOG_LEVEL log_level = (LOG_LEVEL) (verbosity - (int) '0');
                SET_VERBOSITY_LEVEL(log_level);
              }
            }
            break;
          default:
            PrintUsage();
            return -1;
        }
      } else {
        // At this point, we are done processing optional flags, so we
        // get the words.
        word_a = argv[i];
        if( argv[i + 1 ] ) {
          word_b = argv[i + 1];
          i = argc;    // force argument parsing termination
        } else {
          PrintUsage();
          return -1;
        }
      }
      ++i;
    }
  }

  cout << word_a << " " << word_b << endl;

  // This checks that we, indeed, have words to work with in both variables.
  if (!word_a.length() || !word_b.length()) {
    PrintUsage();
    return -1;
  }

  // This rtrims and lowercases the word
  CleanString(word_a);
  CleanString(word_b);

  // We'll now perform our suite of operations on the strings and report results.
  VERBOSE_LOG(LOG_NONE, "Levenshtein:" << endl);
  VERBOSE_LOG(LOG_NONE, CalcLevenshtein(word_a.c_str(), word_b.c_str()) << endl);

  VERBOSE_LOG(LOG_NONE, "Damerau-Levenshtein:" << endl);
  VERBOSE_LOG(LOG_NONE, CalcDamerauLevenshtein(word_a.c_str(), word_b.c_str()) << endl);

  VERBOSE_LOG(LOG_NONE, "Jaro similarity:" << endl);
  VERBOSE_LOG(LOG_NONE, CalcJaroSimilarity(word_a.c_str(), word_b.c_str()) << endl);

  VERBOSE_LOG(LOG_NONE, "Jaro-Winkler similarity:" << endl);
  VERBOSE_LOG(LOG_NONE, CalcJaroWinklerSimilarity(word_a.c_str(), word_b.c_str()) << endl);

  cout << word_a << " " << word_b << endl;

  return 0;
}
