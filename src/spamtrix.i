/* THIS IS A SWIG INTERFACE FILE FOR CREATING WRAPPERS FOR USE */
/* WITH  NON-C++ LANGUAGES*/

%module SpaMtrix
%{
#include "../include/vector.h"
%}
typedef unsigned int idx;
typedef double real;
%include "../include/vector.h"