### He3-simulator
Code for simulating the order parameter of superfluid Helium 3.

Cloned from onsim (Asier Lopez-Eiguren).

Uses HILA framework (Kari Rummukainen et al): https://bitbucket.org/Kari_Rummukainen/hila/src/master/

Create branch materialp for material parameter (Quang Zhang, timohyva@github)

### Compling and ploting the gaps, free energies and bulk phase diagram 
~~~ shellscript
#!/bin/bash
#
...
#

g++ -std=c++17 -g -o pd.app phases_diagram.cpp matep.cpp matep.hpp && ./pd.app

python3 ./pd/sfdata.py
~~~
~~~bash
./sfdata_plots.sh
~~~
