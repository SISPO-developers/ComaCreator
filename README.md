# ComaCreator for Sispo


Dust and gas environment generator for https://github.com/SISPO-developers/sispo. Standalone version has not been tested in a long time so it is mandatory to use Sispo. The plugin is activated by adding keywords (--with-plugin and plugins)

```
   "options": ["--with-render","--with-sim", "--with-plugin"],
    "plugins": ["/app/ComaCreator/interface.py"],
```

where the file path should lead to the root folder of the ComaCreator.

It is recommended to use the https://github.com/SISPO-developers/sispo_docker because that is the most easiest way to use Sispo and it has environment ready for the ComaCreator also.

**How it works**

ComaCreator first generates a map of gas outflow strengths at the surface of the comet. Then the gas outflows are solved around the comet and this data is then used in the particle simulation that takes into account gravity, gas outflows, solar radiation and coriolis forces (if rotating frame). The gas and dust density data is then converted to EXR file that can be loaded by Blender and used for the volumetric scattering density.


##### Requirements

- cython
- pyembree
- trimesh
- tricubic
- numpy
- opencv?

Maybe something else but I can not remember, because so many stuff are already always installed in the Sispo side. 

##### Building and installing

In ""/ComaCreator"

```
python setup.py install
```

##### Running

If the input file is set correctly, the Sispo will invoke JetCreator automatically before rendering (Only Cycles support).

##### Parameters

Inputs are controlled by "defaults.in" which can be created with "create_input.py". Input file in the $ROOT/data/defaults.in overrides the default values of the "$ROOT/src/PROGRAM/root.in". Sispo's input file must contain

```
   "options": ["--with-plugin"],
```

Tell Sispo to use plugins

```
"plugins": ["/app/ComaCreator/interface.py"],
```

Load plugin from /app/ComaCreator/interface.py. This can be also used to call other plugins, but currently no other plugins exits.

```
	"ComaCreator" : {
		"decimate_mesh" : true,
        "decimate_percentage" : 1.0
    }
```

We allow ComaCreator to decimate the shape model (also triangulate) with some percentage (0-1), because it will take forever to iterate 3 million triangle, and it is better to decimate the shape down to 3000 triangles (0.001). 

```
  "comet": {
    "mesh_file": "ROOTPATH/data/comet3.obj"
  },
```

Path to the shape model. Overwritten by the interface. No actions needed.

```
"final_output": "ROOTPATH/data/out.json",
```

Path to the output file that is used to initialize Cycles for Sispo. No actions needed

```
  "JetCreator": {
    "ID": "myID",
    "programID": {
      "NSDSGenerator": "src.NSDSGenerator.wrapper",
      "KGasSimulation": "src.KGasSimulation.wrapper",
      "KParticleSimulation": "src.KParticleSimulation.wrapper",
      "ToBlender": "src.ToBlender.wrapper"
    },
    "common": {
      "executionID": "THIS_IS_HOLDER_DO_NOT_USE",
      "write_output_to": "default",
      "output_folder": "ROOTPATH/data/",
      "nProcesses": 3,
      "tmpFOLDER": "ROOTPATH/tmp/",
      "special_instruction": "empty"
    },
```

"nProcesses", controls the number of invoked processes. So if you have 16 cores in your CPU, you might want to set this to 15.

Otherwise, no actions required. These are used to control the program flow. ProgramID is used to specify what kidn of programs the pipeline consists. So if someone wants to expand the pipeline it is possible to add "TESTADDITION" : "src.testaddition.wrapper" and then the interface.py would try to find input parameters for "TESTADDITION" from the existing input file (see below) and execute src.testaddition.wrapper. 

    "sim_common": {
      "mesh_file": "default",
      "T_gas": 210.0,
      "m_gas": 0.01801528,
      "R_gas": 0.00046152
    },
T_gas, is the temperature of the gas, m_gas and R_gas is the mass and the gas constant

    "pipeline": [
      "DefaultGasPreprocessor",
      "DefaultGasDensityCreator",
      "DefaultParticleGenerator",
      "DefaultToBlender"
    ],
Specify the pipeline. The programs are run in this order

    "DefaultGasPreprocessor": {
      "programID": "NSDSGenerator"
    },
Inputs for preprocessor. See (src/PROGRAM/defaults.in) for more input parameters. 

    "DefaultGasDensityCreator": {
      "programID": "KGasSimulation",
      "outflow_rates": "prevOUTPUT",
      "nBounds": 3,
      "mesh_scale" : 1.0,
      "sample_points": [
        64,
        64,
        64
      ],
      "bbox_lower": [
        -5.0,
        -5.0,
        -5.0
      ],
      "bbox_upper": [
        5.0,
        5.0,
        5.0
      ],
      "include_data": [
        "sim_common"
      ],
      "special_instruction": "empty"
    },
Inputs for the gas environment generator. The most important stuff is the sample_points. This is the resolution of the space around the comet. bbox_lower, and bbox_upper are the bounding box limits. Use cubes ([-x,-x,-x], [x,x,x]), allowing more varied shape was a legacy option. See (src/PROGRAM/defaults.in) for more input parameters. 

    "DefaultParticleGenerator": {
      "mesh_scale" : 1.0,
      "programID": "KParticleSimulation",
      "gas_field": "prevOUTPUT",
      "include_data": [
        "sim_common"
      ],
      "n_particles": 10000,
      "special_instruction": "empty"
    },
Inputs for the dust environment generator. "n_particles" controls how many particles are simulated. See (src/PROGRAM/defaults.in) for more input parameters. 

    "DefaultToBlender": {
      "programID": "ToBlender",
      "include_data": [
        "sim_common"
      ],
      "write_output_to": "final",
      "dust_particle_field": "prevOUTPUT"
    }
Conversion from ComaCreator to EXR. See (src/PROGRAM/defaults.in) for more input parameters

##### Based on 

- Kramer, Tobias & Noack, Matthias & Baum, Daniel & Hege, Hans-Christian & Heller, Eric. (2018). Dust and gas emission from cometary nuclei: the case of comet 67P/Churyumov–Gerasimenko. Advances in Physics: X. 3. 1404436. 10.1080/23746149.2017.1404436. 
- Tobias Kramer, Matthias Läuter, Martin Rubin, Kathrin Altwegg, Seasonal changes of the volatile density in the coma and on the surface of comet 67P/Churyumov–Gerasimenko, Monthly Notices of the Royal Astronomical Society, Volume 469, Issue Suppl_2, July 2017, Pages S20–S28, https://doi.org/10.1093/mnras/stx866

But the ComaCreator doesn't function as they intended, has more bugs certainly, and we have no collaboration with them, so please don't use ComaCreator as substitute for a real science. ComaCreator is in development stage maybe forever and should be only used for creating pretty realistic looking pictures (about the realistic and pretty we can argue of course :D). 

**Want to help**

* How to increase performance. 

* Increase resolution (nowadays OpenVDB is available so could we use this)

* Reduce file size (We use EXR and the file size can grow to 1 gb easy)