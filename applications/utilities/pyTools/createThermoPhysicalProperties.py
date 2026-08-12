import sys


species = sys.argv[1].split(",")
diffusion = sys.argv[2].split(",")
surface_species = sys.argv[3].split(",")
surface_masters = sys.argv[4].split(",")
eps_chemistry_min = sys.argv[5] if len(sys.argv) > 5 else "0.0"

if len(species) != len(diffusion):
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
    f"epsChemistryMin {eps_chemistry_min};\n",
    "kineticPhases\n",
    "{\n",
    "}\n",
    "kineticPhaseReactions\n",
    "{\n",
    "}\n",
    "solutionSpecies\n",
    "{\n",
]

for name, coefficient in zip(species, diffusion):
    lines.extend(
        [
            f"    {name}\n",
            "    {\n",
            f"        D D [ 0 2 -1 0 0 0 0 ] {coefficient};\n",
            "    }\n",
        ]
    )

lines.extend(["};\n", "\n", "surfaceSpecies\n", "{\n"])

for name in surface_species:
    lines.append(f"    {name};\n")

lines.extend(["};\n", "\n", "surfaceMasters\n", "{\n"])

for name in surface_masters:
    lines.append(f"    {name};\n")

lines.extend(["};\n", "\n", "// ************************************************************************* //\n"])

with open("constant/thermoPhysicalProperties", "w", encoding="utf-8") as output:
    output.writelines(lines)
