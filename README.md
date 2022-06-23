# Satellite Network Simulator
This code is used as an end-to-end satellite network simulation environment designed for verification purposes. It is meant to be an easy-to-use open-source simulation tool for anyone dealing with satellite communication services, or just willing to contribute with software development.

## Files
The code is based in object-oriented programming and it is structured in three different scripts, i.e. <i>main.py</i>, <i>simulator.py</i> and <i>objects.py</i>. Each of the scripts have different levels of functionality. The simulation is executed via command line using the <i>main.py</i> script, where the default parameters of the simulation can be modified using argument parsing. Call `python3.5 ./main.py -h` to display options.

### Arguments
|   Flag   |                                   Description                                                  |
|----------|:----------------------------------------------------------------------------------------------:|
|   `-h`   |   display argument options             |
|   `-v`   |   [y] enables 3d visualization         |
|   `-dt`  |   simulation sampling time in hours    |
|   `-c`   |   followed by a .json file with the parameters of the satellite constellation    |
|   `-gs`  |   followed by a .json file with the parameters of the ground stations   |
|   `-sh`  |   followed by a .json file with the parameters of the ships    |
|   `-a`   |   followed by a .json file with the parameters of the sea areas   |

E.g.: 
`python3.5 ./main.py -v y -dt 0.001 -c constellation.json -gs gs.json -sh ships.json -a areas.json`

## Mission Customization
The customization of the mission can be performed in <i>simulator.py</i> by calling the desired high-level functionality inside the simulation loop. You can find the low-level code and the available functions for each object in <i>objects.py</i>.


## Source Code
### Preparation and requirements
The simulator consists of the following 3 scripts:
<ol>
  <li><i>main.py</i>
  <li><i>simulator.py</i> 
  <li><i>objects.py</i>
</ol>

and the following 4 input json files:
<ol>
  <li><i>constellation.json</i>
  <li><i>gs.json</i>
  <li><i>ships.json</i>
  <li><i>areas.json</i>
</ol>

The libraries required are:
<ul>
  <li><i>math</i>
  <li><i>numpy</i> 
  <li><i>scipy</i>
  <li><i>pyqtgraph</i>
</ul>

It might be that extra libraries are needed to import pyqtgraph. The following are some that might be needed, if not installed yet.
<ul>
  <li><i>libxcb-randr0-dev</i>
  <li><i>libxcb-xtest0-dev</i> 
  <li><i>libxcb-xinerama0-dev</i>
  <li><i>libxcb-shape0-dev</i>
  <li><i>libxcb-xkb-dev</i>
</ul>
