#define _USE_MATH_DEFINES
#define USE_ASCENT 
#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
//#include <math.h>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"
#include "pario.hpp"

//#if defined USE_ASCENT
#include "ascent.hpp"
#include "conduit_blueprint.hpp"
//#endif

void parIO::shutdown() {
    pio.close();
}
