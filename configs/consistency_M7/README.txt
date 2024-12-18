Configuration files to use for validating cell densities of the M7 circuit

Note: In all files the paths to the input files, i.e. the density .nrrd files, must still be updated to their actual locations!
Similarly, for the validation of densities in the circuit, the path to the circuit needs to be updated!

Validation of consistency of "input" densities
CONFIG: consistency_checks_density_inputs.json

This refers to the densities generated by Dimitri's algorithms, i.e. VIP, SST, PV, etc.; but also various glia types.
I compiled a list of various logical consistency checks.
I call them "input" densities, because these files are the input for the next step, generation of "output" densities.

-----------

Validation of consistency of inhibitory "output" densities
CONFIG: consistency_checks_inh_density_outputs.json

This refers to the densities generated by Yann's algorithms, i.e. individual inhibitory mtypes, such as L4_MC.
I compiled a list of logical consistency checks, such that inhbitory types add up to the prescribed "GAD" density.
I added separate checks for each cortical region, i.e. we validate that the constraints are not just fulfilled on average,
but also in subregions.

------------

Validation of consistency of excitatory "output" densities
CONFIG: consistency_checks_exc_density_outputs.json

This refers to the densities generated by Hugo's algorithms, i.e. individual excitatory mtypes, such as L4_MC.
It's essentially the same as above, just for excitatory instead of inhibitory m-types. Validation also region-per-region.

------------

Validation of placed cells
CONFIG: circuit_test_all_densities.json

This validates that the densities of actually placed neurons of various types match the prescriptions.
Validations are performed region-per-region and for individual m-types. But not both at the same time, as that would lead
to a rather large number of individual tests.
