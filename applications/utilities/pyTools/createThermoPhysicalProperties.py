import sys

species=sys.argv[1].split(',')
Diff=sys.argv[2].split(',')
SurfSpecies=sys.argv[3].split(',')
SurfMasters=sys.argv[4].split(',')

if len(species) != len(Diff):
    raise ValueError(
        "The species and diffusion coefficient tables must have the same length"
    )

lines = [
    "FoamFile\n",
    "{\n",
    "    version     2.0;\n",
    "    format      ascii;\n",
    "    class       dictionary;\n",
    '    location    "constant";\n',
    "    object      thermoPhysicalProperties;\n",
    "}\n",
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
    "\n",
    "mixtureType phreeqcMixture;\n",
    "kineticPhases\n",
    "{\n",
    "}\n",
    "kineticPhaseReactions\n"
    "{\n"
    "}\n"

    "solutionSpecies\n",
    "{\n",
]

for i in range(len(species)):
    lines.extend(
        [
            f"    {species[i]}\n",
            "    {\n",
            f"        D D [ 0 2 -1 0 0 0 0 ] {Diff[i]};\n",
            "    }\n",
        ]
    )

lines.extend(["};\n", "\n", "surfaceSpecies\n", "{\n"])

for i in range(len(SurfSpecies)):
    lines.append(f"    {SurfSpecies[i]};\n")

lines.extend(["};\n", "\n", "surfaceMasters\n", "{\n"])

for i in range(len(SurfMasters)):
    lines.append(f"    {SurfMasters[i]};\n")

lines.extend(["};\n", "\n", "// ************************************************************************* //\n"])

with open("constant/thermoPhysicalProperties", "w") as file:
    file.writelines(lines)
