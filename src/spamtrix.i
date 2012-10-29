/* THIS IS A SWIG INTERFACE FILE FOR CREATING WRAPPERS FOR USE */
/* WITH  NON-C++ LANGUAGES*/

%module SpaMtrix
%{
#include "../include/setup.h"
#include "../include/vector.h"
#include "../include/ircmatrix.h"
#include "../include/matrixmaker.h"
#include "../include/fleximatrix.h"
%}
typedef unsigned int idx;
typedef double real;
%include "../include/setup.h"
%include "../include/vector.h"
%include "../include/ircmatrix.h"
%include "../include/matrixmaker.h"
%include "../include/fleximatrix.h"