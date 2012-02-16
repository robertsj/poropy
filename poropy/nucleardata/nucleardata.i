/* nucleardata.i */
%module nucleardata
%{
/* Includes the header in the wrapper code */
#include "nucleardata.h"
#include "nucleardata_IFBA.h"
#include "nucleardata_WABA.h"
#include "nucleardata_GAD.h"
%}

/* Parse the header file to generate wrappers */
%include "nucleardata.h"
%include "nucleardata_IFBA.h"
%include "nucleardata_WABA.h"
%include "nucleardata_GAD.h"
