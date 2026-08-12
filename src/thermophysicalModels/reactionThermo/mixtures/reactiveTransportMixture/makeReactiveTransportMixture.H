/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

InClass
    Foam::reactiveTransportMixture

Description
    Macro for registering reactiveTransportMixture template instantiations.

\*---------------------------------------------------------------------------*/

#ifndef makeReactiveTransportMixture_H
#define makeReactiveTransportMixture_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactiveTransportMixture(MixtureType)                             \
                                                                              \
typedef reactiveTransportMixture<MixtureType>                                 \
    reactiveTransportMixture##MixtureType;                                    \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    reactiveTransportMixture##MixtureType,                                    \
    "reactiveTransportMixture<" #MixtureType ">",                             \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicReactiveTransportMixture,                                            \
    reactiveTransportMixture##MixtureType,                                    \
    fvMesh                                                                    \
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
