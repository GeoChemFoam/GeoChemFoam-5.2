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

\*---------------------------------------------------------------------------*/

#include "basicReactiveTransportMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicReactiveTransportMixture, 0);
    defineRunTimeSelectionTable(basicReactiveTransportMixture, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicReactiveTransportMixture::basicReactiveTransportMixture
(
    const fvMesh& mesh
)
:
    basicMultiComponentTransportMixture(mesh),
    kineticPhase_(),
    kineticPhaseReaction_(),
    Mws_(),
    rhos_(),
    Ke_(),
    k0_(),
    scoeff_(),
    kiSpeciesIndex_(),
    ki_(),
    Ws_(),
    Omega_(),
    R_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::basicReactiveTransportMixture::~basicReactiveTransportMixture()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicReactiveTransportMixture>
Foam::basicReactiveTransportMixture::New(const fvMesh& mesh)
{
    word mixtureTypeName;

    // Keep this temporary out of the objectRegistry before the selected
    // mixture constructs its own thermoPhysicalProperties dictionary.
    {
        IOdictionary thermoDict
        (
            IOobject
            (
                "thermoPhysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        thermoDict.lookup("mixtureType") >> mixtureTypeName;
    }

    const word reactiveTransportMixtureType
    (
        "reactiveTransportMixture<" + mixtureTypeName + ">"
    );

    Info<< "Selecting reactiveTransportMixture type "
        << reactiveTransportMixtureType << endl;

    auto cstrIter =
        fvMeshConstructorTablePtr_->find(reactiveTransportMixtureType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reactiveTransportMixture type "
            << reactiveTransportMixtureType << nl << nl
            << "Valid reactiveTransportMixture types are:" << nl
            << fvMeshConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<basicReactiveTransportMixture>(cstrIter()(mesh));
}


// ************************************************************************* //
