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

#include "reactiveTransportMixture.H"
#include "inertMultiComponentMixture.H"
#include "phreeqcMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::reactiveTransportMixture<MixtureType>::reactiveTransportMixture
(
    const fvMesh& mesh
)
:
    basicMultiComponentTransportMixture(mesh),
    basicReactiveTransportMixture(mesh),
    MixtureType(*this, this->subDict("solutionSpecies").toc(), mesh),
    component1Index_(),
    component2Index_()
{
    kineticPhase_ = this->subDict("kineticPhases").toc();
    kineticPhaseReaction_ = this->subDict("kineticPhaseReactions").toc();

    Mws_.setSize(kineticPhaseReaction_.size());
    rhos_.setSize(kineticPhaseReaction_.size());
    Ke_.setSize(kineticPhase_.size());
    k0_.setSize(kineticPhaseReaction_.size());
    scoeff_.setSize(kineticPhaseReaction_.size());
    kiSpeciesIndex_.setSize(kineticPhaseReaction_.size());
    ki_.setSize(kineticPhaseReaction_.size());
    Ws_.setSize(kineticPhase_.size());
    Omega_.setSize(kineticPhase_.size());
    R_.setSize(kineticPhaseReaction_.size());


    component1Index_.assign(kineticPhase_.size(), -1);
    component2Index_.assign(kineticPhase_.size(), -1);

    const dictionary& kpDict = this->subDict("kineticPhases");

    forAll(kineticPhase_, i)
    {
        const dictionary& kpSubDict = kpDict.subDict(kineticPhase_[i]);
        speciesTable kpSpecies(kpSubDict.subDict("species").toc());

        forAll(this->species_, j)
        {
            if (this->species_[j] == kpSpecies[0])
            {
                component1Index_[i] = j;
                break;
            }
            if (j == this->species_.size())
            {
                Info<< kpSpecies[0] << "not find in solutionSpecies"
                    << endl
                    << abort(FatalError);
            }
        }

        forAll(this->species_, j)
        {
            if (this->species_[j] == kpSpecies[1])
            {
                component2Index_[i] = j;
                break;
            }
            if (j == this->species_.size())
            {
                Info<< kpSpecies[1] << "not find in solutionSpecies"
                    << endl
                    << abort(FatalError);
            }
        }

        Ke_.set
        (
            i,
            new dimensionedScalar("Ke", kpSubDict)
        );

        Ws_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Ws_" + kineticPhase_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Ws_" + kineticPhase_[i], dimless, 1)
            )
        );

        Omega_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    kineticPhase_[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->Y_[component1Index_[i]]*this->Y_[component2Index_[i]]/Ke_[i]
            )
        );
    }

    const dictionary& kprDict = this->subDict("kineticPhaseReactions");

    forAll(kineticPhaseReaction_, j)
    {
        const dictionary& kprSubDict = kprDict.subDict(kineticPhaseReaction_[j]);

        Mws_.set
        (
            j,
            new dimensionedScalar("Mws", kprSubDict)
        );

        rhos_.set
        (
            j,
            new dimensionedScalar("rhos", kprSubDict)
        );

        k0_.set
        (
            j,
            new dimensionedScalar("k0", kprSubDict)
        );

        scoeff_.set
        (
            j,
            new scalarList(this->species_.size(), 0.0)
        );

        const dictionary& kprSpeciesDict = kprSubDict.subDict("species");
        speciesTable kprSpecies(kprSpeciesDict.toc());

        forAll(this->species_, i)
        {
            if (kprSpeciesDict.isDict(this->species_[i]))
            {
                scoeff_[j][i] =
                    readScalar(kprSpeciesDict.subDict(this->species_[i]).lookup("scoeff"));
            }
        }

        kiSpeciesIndex_.set
        (
            j,
            new labelList(kprSpecies.size(), 0)
        );

        forAll(kprSpecies, i)
        {
            forAll(this->species_, k)
            {
                if (kprSpecies[i] == this->species_[k])
                {
                    kiSpeciesIndex_[j][i] = k;
                }
            }
        }

        ki_.set
        (
            j,
            new PtrList<dimensionedScalar>(kprSpecies.size())
        );

        forAll(kprSpecies, i)
        {
            ki_[j].set
            (
                i,
                new dimensionedScalar("ki", kprSpeciesDict.subDict(kprSpecies[i]))
            );
        }

        R_.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "R_" + kineticPhaseReaction_[j],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                k0_[j]
            )
        );
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::reactiveTransportMixture<MixtureType>::~reactiveTransportMixture()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::reactiveTransportMixture<MixtureType>::correct()
{
    MixtureType::correct();

    forAll(kineticPhase_, i)
    {
        Omega_[i] =
            this->Y_[component1Index_[i]]*this->Y_[component2Index_[i]]/Ke_[i];
    }

    forAll(kineticPhaseReaction_, j)
    {
        R_[j] = k0_[j];

        forAll(kiSpeciesIndex_[j], i)
        {
            R_[j] += ki_[j][i]*this->Y_[kiSpeciesIndex_[j][i]];
        }
    }
}


// ************************************************************************* //

namespace Foam
{
    template class reactiveTransportMixture<inertMultiComponentMixture>;
    template class reactiveTransportMixture<phreeqcMixture>;
}
