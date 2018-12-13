#ifndef _ANAGRAM_FLAGS_H_
#define _ANAGRAM_FLAGS_H_

// StringbenchFlags
// This structure defines a 32-bit word of flags
union StringBenchFlags {
  struct {
    unsigned int tree_engine : 1;
    unsigned int allow_dupes : 1;
    unsigned int output_directly : 1;
    unsigned int big_dictionary: 1;
    unsigned int print_subset : 1;
  } bit;
  unsigned int word;
};

#endif // #ifndef _ANAGRAM_FLAGS_H_
