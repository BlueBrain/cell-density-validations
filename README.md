# Consistency test specifications
The configuration file is a json file specifying a dict with the following entries:
## Required entries
  - annotations: Value is a path to a .nrrd file. The file holds the brain region annotation volume
  - hierarchy: Value is a path to a .json file. The file holds the brain region hierarchy and ids.
  - nrrd: Value is a dict. Each value of the dict is a path to an .nrrd file with cell densities, 
    or a list of such files. Several files can be grouped together in the same list if they are used together in all
    tests specified below. The keys of the dict are labels given to the files (nrrd-labels) to be used in the rest 
    of the specifications
  - tests: Value is a list of tests, formatted as below

## Optional entries
  - circuits: Value is a dict. Each value of the dict is the path to a CircuitConfig file that can be used with Bluepy.
The keys are used as shorthand names for the circuits in "circuit-target" specs (see below). If you only want to validate atlases,
    and not circuits, then you can leave this out.

### Specification of individual tests
Individual tests are specified as dicts with the following entries:
#### Required entries
  - left: A list of nrrd-labels or numbers or "circuit-target" specs (see below)
  - right: Another list of nrrd-labels or numbers or "circuit-target" specs (see below)
The lists are specified such that sum of the volumes in "left" is expected to be (approximately) equal to the sum of volumes in "right".
If a number is part of any of the lists, it will be added as a constant.
#### Optional entries
  - tolerance: A dict with any number of the following keys: ["absolute", "relative"].
    Values specify the maximal error tolerated in terms of absolute or relative difference, as given by the key. 
    If both "absolute" and "relative" are used, then both tests must pass.
    If not used, then no tolerance is given, i.e. "left" must exactly equal "right"
  - mask: A specification of brain regions. The test will be applied only to regions matching the description.
    If no mask is given, all voxels are used. Specification of regions: See below.
The test will always compare the total sums of voxels on the "left" and "right" and check whether their difference exceeds tolerance.
If only nrrd-labels are used in both "left" and "right", additionally a per-voxel comparison of the sum on the left and right are performed.
The per-voxel test checks whether the maximum error over (masked) voxels exceeds tolerance.

#### Mask specification
A mask specification one of:
  - A string that is the acronym of a brain region used in the specified hierarchy.json file (see above)
  - A string starting with "Layer" and the rest of it specifies a cortical layer using an arabic numeral.
    Example: "Layer5". A range of layers is specified as follows: "Layer[2-4]"
  - A list of length 2:
    - First entry is one of "and", "or"
    - Second entry is a list of valid mask specifications. 
  If "and" is used, then a voxel has to belong to all specifications in the list, 
  if "or" is used, then a voxel only has to belong to at least one.

#### Circuit-target specs
This specifies a set of neurons in a circuit model, where the properties of the neurons fulfill certain, specifiable criteria.
When used as part of an individual test (see above), the density of these neurons in the context of the annotation atlas (see above under "required entries") is used for validation.
That means, a circuit-target spec can be used in exactly the same way as an .nrrd file for the purpose of validations, but it will use the density of actually placed neurons!

A "circuit-target" spec is a dict with the following entries:
  - "circuit": The value must be a string that can be used as index of the "circuits" dict specified above under "optional entries".
  - "group" (Optional): When provided, this limits the neurons considered to the specified cell target. The target must be defined in the specified circuit.
  - "filters" (Optional): An optional dict. Keys are names of neuron properties, each value is either a list of valid values or simply a single valid value.
Then, only neurons that match for all given properties the valid value(s) are considered. Potential properties are the ones that can be used in 
    bluepy when calling Circuit.cells.get(properties=[...]). The exact list depends on the circuit, the most common are:
    - 'etype'
    - 'layer'
    - 'me_combo' 
    - 'morph_class'
    - 'morphology'
    - 'mtype'
    - 'region'
    - 'synapse_class'
    
Example: {"Circuit": "my_circ", "group": "Excitatory"} and {"Circuit": "my_circ", "filters": {"synapse_class": "EXC"}}
should both result in the same group of neurons.