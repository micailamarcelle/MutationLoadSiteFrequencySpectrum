#ifndef GLOBAL_VARS_H_INCLUDED
#define GLOBAL_VARS_H_INCLUDED 1

#define VERBOSE 0
#define VERYVERBOSE 0
#define MISCELLANEOUS 1
#define INDIVIDUALWIDATA 1
#define OVERLAPPINGGENERATIONS 1 //probably should be an input argument at some point.
#define LSB(i) ((i) & -(i)) //isolates least significant single bit for fenwick tree
#define PI 3.141592654

#define RENORM_FREQ 500 // specifies how often (in terms of the number of generations) renormalization should occur
#define CAPTURE_FREQ 5000 // specifies how often (in terms of the number of generations) the tskit data should be captured

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }//error checking for tree sequence recording

#endif // GLOBAL_VARS_H_INCLUDED

