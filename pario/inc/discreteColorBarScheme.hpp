#ifndef DISCRETECOLORBARSCHEME_H
#define DISCRETECOLORBARSCHEME_H

#include <vector>
#include "conduit_blueprint.hpp"

// Jet 9 colors
float colors[9][3] = {
  {0.0, 0.0, 0.5}, //pM 1
  {0.0, 0.0, 1.0}, //pM 2
  {0.0, 0.5, 1.0}, //pM 3
  {0.0, 1.0, 1.0}, //pM 4
  {0.5, 1.0, 0.5}, //pM 5
  {1.0, 1.0, 0.0}, //pM 6
  {1.0, 0.5, 0.0}, //pM 7
  {1.0, 0.0, 0.0}, //pM 8
  {0.5, 0.0, 0.0}  //pM 9
};

void colorControlPointsGenerator(conduit::Node &control_points) {
 for(unsigned int i = 0; i < 9; i++)
  {
    double p0 = double(i) / 9.0;
    double p1 = double(i+1) / 9.0;

    // left edge
    {
        conduit::Node &cp = control_points.append();
        cp["type"] = "rgb";
        cp["position"] = p0;
        cp["color"].set(std::vector<float>{
            colors[i][0],
            colors[i][1],
            colors[i][2]
        });
    }

    // right edge (same color!)
    {
        conduit::Node &cp = control_points.append();
        cp["type"] = "rgb";
        cp["position"] = p1 - 1e-6;   // - 1e-6 handle the precession
        cp["color"].set(std::vector<float>{
            colors[i][0],
            colors[i][1],
            colors[i][2]
        });
    }
  } // for loop ends here

} // colorControlPointGenerator() ends here
#endif
