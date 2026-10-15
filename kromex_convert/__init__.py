"""kromex_convert: external data as KROME-X inputs, part of the KROME-X
install (independent of the code generation: they write files it reads):

    kida2krome   a KIDA network (kida.uva)        -> .knet
    umist2krome  a UMIST ratefile (RATE22, RATE12) -> .knet
    lamda2krome  a LAMDA molecular data file       -> coolant file
    chianti2krome  CHIANTI ions                    -> coolant file

Each is a command (installed with pip) and, from a checkout, python3 -m
kromex_convert.NAME. docs/manual.md section 2 and the [thermal] options."""
