# SPDX-License-Identifier: Apache-2.0
import voxcell
import os

def find_file_at(fn, root):
    if not os.path.isfile(fn):
        fn = os.path.join(root, fn)
    return fn

def load_with_arithmetics(specs, at_root=None):
    if isinstance(specs, tuple) or isinstance(specs, list):
        operation = specs[0]
        operands = specs[1]
        operand_vals = [load_with_arithmetics(_x, at_root=at_root) for _x in operands]

        if operation == "divide":
            assert len(operand_vals) == 2
            return operand_vals[0] / operand_vals[1]
        if operation == "subtract":
            assert len(operand_vals) == 2
            return operand_vals[0] - operand_vals[1]
        if operation == "add":
            v = 0
            for oprnd in operand_vals:
                v = v + oprnd
    
    if isinstance(specs, str):
        fn = find_file_at(specs, at_root)
        v = voxcell.VoxelData.load_nrrd(fn)
        return v.raw
    return specs
