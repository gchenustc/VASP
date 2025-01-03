# !/bin/env bash
# Desc: Make super cell for $1 parameter, For Example, ./script.sh POSCAR_name 3 3 3 means expand the POSCAR_name for (3,3,3) super cell. The Multiwfn software is must.
# Time: 2024-12-29
Multiwfn $1 << EOF >/dev/null
300
7
19
$2
$3
$4
-10
0
100
2
27
POSCAR_$2-$3-$4
0
q
EOF
