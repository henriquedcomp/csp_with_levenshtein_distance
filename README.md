# Closest String Problem with Levenshtein Distance

* [1. Installing OR-Tools](#1-installing-or-tools)
* [2. Overview and project structure](#2-overview-and-project-structure)
* [3. Instances](#3-instances)
* [4. How to run](#4-how-to-run)

<hr>

## 1. Installing OR-Tools
In order to run this project, you'll need OR-Tools v9.14 with Gurobi v12.0.3. This combination certainly will work together [\[Gurobi not working in v9.15\]](https://github.com/google/or-tools/discussions/5026#discussioncomment-15708012). Additionally, C++ Tools is necessary.

In this section I'll help you to install OR-Tools. Gurobi can easily be setup following [this article](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer).

Before we proceed, I'll suppose, once you've arrived at this repository, you're a Linux user. If it's not the case, you should consider switch to Linux or follow documentation for your OS. [OR-Tools for C++](https://developers.google.com/optimization/install/cpp).

First, let's clone the OR-Tools repository. I cloned it into my HOME directory.
```
git clone https://github.com/google/or-tools.git
cd or-tools
```

Inside the repository, switch to release 9.14
```
git checkout v9.14
```

Next, let's generate the OR-Tools build.
```
cmake -S . -B build -DBUILD_DEPS=ON
```

The next step will take some time. We are going to compile the source code. The -j$(nproc) flag will locate CPU cores on your machine and use then to perform parallel processing. In my case, this led to a crash, so I have set a limit of 6 cores.
Command using all cores:
```
cmake --build build --config Release -j$(nproc)
```

Command limiting to 6 cores:
```
cmake --build build --config Release -j6
```

It can print some warnings, this is normal. It took about 30 to 40min in my case.

Finally, we'll install the OR-Tools. By default, it will be located in /usr/local
```
sudo cmake --install build
```

IMPORTANT: Before running this project, if you are using a python environment, please exit it. Such environments may change some important environment variables and cause conflicts when running OR-Tools.
For example, if you're using Miniconda, exit the environment by using:
```
conda deactivate
```

Now, OR-Tools must be working properly with Gurobi in this C++ project.

## 2. Overview and project structure
This project aims to test some integer linear programming formulations to closest string problem (aka center string problem) using Levenhtein distance as a distance metric. Also, we'll improve formulations by adding some new constraints and make tests using constraint programming (CP-SAT solver) too.

The project is organized as follows:
```
csp_with_levenshtein_distance/
├── include/               # Folder with .hpp reusable files
│   ├── config.hpp         # File for code global configs
│   └── {...}
├── instances/             # Problem instances
├── results/               # Instance execution results
├── src/                   # Main project's source code
│   ├── main.cpp           # Program's entry point
│   └── {...}             
├── .gitignore
├── CMakeLists.txt         # CMake Project file
├── generate_instances.sh
├── instance_generator.cpp
├── README.md
├── run.sh
└── runall.sh
```

You can change execution time limit and Levenshtein distance operation costs in include/config.hpp.
Attention: In order to keep compatibility with constraint programming, use only integer costs.
## 3. Instances
Random instances were created in such a manner that you can easily rebuilt it with other parameters as needed. The file instance_generator.cpp owns this responsability.

Before run generate_instances script (that calls instance_generator), make sure it has execution permission:
```
chmod +x generate_instances.sh
```

## 4. How to run
On the first run:
```
cmake -S . -B build
```
This will create the build folder.

Make sure also that shellscript files had execution permission:
```
chmod +x run.sh
chmod +x runall.sh
```
Now, you can just run _run.sh_ to compile code and execute an instance or _runall.sh_ to compile code and run all instances. You can select an instance to run in _run.sh_ writing its path in the script.
