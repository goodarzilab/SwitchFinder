import numpy as np

_stem = np.uint8(1)
_loop = np.uint8(2)

_left_stem = np.uint8(3)
_right_stem = np.uint8(4)

_T = np.uint8(1)
_C = np.uint8(2)
_G = np.uint8(3)
_A = np.uint8(4)
_none = np.uint8(5)

nt_list = [_T, _C, _G, _A]

_base_pairing_dict = {_T : _A,
                      _A : _T,
                      _G : _C,
                      _C : _G
}

_char_to_nt_mapping =  {_T : 'T',
                        _C : 'C',
                        _G : 'G',
                        _A : 'A',
                       }


_nt_to_char_mapping = {'T' : _T,
                       'C' : _C,
                       'G' : _G,
                       'A' : _A
                       }


_char_to_numpy = {
                    "(" : _left_stem,
                    "<" : _left_stem,
                    "." : _loop,
                    ")" : _right_stem,
                    ">" : _right_stem
                    }

_char_to_structure = {
                    "(" : _stem,
                    "<" : _stem,
                    "." : _loop,
                    ")" : _stem,
                    ">" : _stem
                    }

STRUCT_LIST = np.array([_stem, _loop])
STRUCT_LIST_CHAR = ['stem', 'loop']


_char_to_extended_structure = {
                    "(" : _left_stem,
                    "<" : _left_stem,
                    "." : _loop,
                    ")" : _right_stem,
                    ">" : _right_stem
                    }

_extended_structure_to_char = {
    _left_stem : "<",
    _right_stem : ">",
    _loop : "."
}

# for making pallettes for SHAPE data
_SHAPE_VMIN = 0
_SHAPE_VMAX = 1
