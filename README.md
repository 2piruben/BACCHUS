# BACCHUS
## **BAC**teria **C**olony **H**ybrid **U**nified **S**imulator

_(Some people will ask "Unified?", Bacchus will answer: Wine not?)_

---

Bacchus is a hybrid agent-based model simulating bacterial colonies in which each bacteria can be provided by an internal set of reactions and bacteria can interact by secreting diffusible molecules. 


# Files


- **`bacterium.cpp`, `bacterium.h`** the bacterium class contains the properties of each single bacterium in the population
-  **`population.cpp`, `population.h`** the population class manages groups of bacterium objects controlling their evolution in time
-  **`algebra2d.cpp`, `algebra2d.h`** library to perform vector operations in the plane
-  **`pars.in`** input parameter file
-  **`main.cpp`** main file to run the simulation


# TODO list:

- Introduce a class to manage diffusible molecules
- Introduce a class to manage intracellular reactions
- Change the ".in" file to a ".xml" or ".json"