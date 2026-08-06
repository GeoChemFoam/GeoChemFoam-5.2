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

#include "basicMultiComponentTransportMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMultiComponentTransportMixture::basicMultiComponentTransportMixture
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "thermoPhysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    DY_(this->subDict("solutionSpecies").toc().size()),
    hasTensorDY_(this->subDict("solutionSpecies").toc().size()),
    DTensorY_(this->subDict("solutionSpecies").toc().size())
{
    wordList specieNames(this->subDict("solutionSpecies").toc());

    Info<< "Read species diffusion coefficients\n" << endl;

    forAll(specieNames, i)
    {
        const dictionary& solutionSpeciesDict = this->subDict("solutionSpecies");
        const dictionary& subdict = solutionSpeciesDict.subDict(specieNames[i]);

        DY_.set
        (
            i,
            new dimensionedScalar("D", subdict)
        );

        const word DTensorName("D_" + specieNames[i]);

        hasTensorDY_[i] = IOobject
        (
            DTensorName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ).typeHeaderOk<volTensorField>();

        if (hasTensorDY_[i])
        {
            Info<< "Reading tensorial DTensor for species " << specieNames[i]
                << " from " << mesh.time().timeName() << "/" << DTensorName
                << nl << endl;

            DTensorY_.set
            (
                i,
                new volTensorField
                (
                    IOobject
                    (
                        DTensorName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            Info<< "No tensorial DTensor found for species " << specieNames[i]
                << " - using scalar DY only" << nl << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::basicMultiComponentTransportMixture::~basicMultiComponentTransportMixture()
{}


// ************************************************************************* //
