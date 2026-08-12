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

#include "phreeqcMixture.H"
#include "RM_interface_C.h"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"
#include "reactingWallFvPatchScalarField.H"
#include "error.H"
#include "Pstream.H"
#include <algorithm>
#include <chrono>
#include <vector>
#include <iomanip>

namespace Foam
{

void phreeqcMixture::checkPhreeqcStatus
(
    const char* operation,
    const int status
) const
{
    if (status >= 0)
    {
        return;
    }

    const int decoded = phreeqcRM_->decodeError(status);
    const std::string backendError = phreeqcRM_->getErrorString();

    Ostream& os = FatalErrorInFunction;
    os
        << "PHREEQC operation '" << operation << "' failed." << nl
        << "status=" << status << nl
        << "decodedStatus=" << decoded;

    if (!backendError.empty())
    {
        os << nl << "backendError:" << nl << backendError;
    }

    os << exit(FatalError);
}

void phreeqcMixture::configureSurfaceMode(const dictionary& thermoDict)
{
    surfaceSolutionIndex_.assign(static_cast<std::size_t>(mesh_.nCells()), -1);

    const bool hasSurfaceSpecies = thermoDict.found("surfaceSpecies");
    const bool hasSurfaceMasters = thermoDict.found("surfaceMasters");
    const bool emptySurfaceSpecies =
        hasSurfaceSpecies && thermoDict.subDict("surfaceSpecies").toc().empty();
    const bool emptySurfaceMasters =
        hasSurfaceMasters && thermoDict.subDict("surfaceMasters").toc().empty();

    if
    (
        !(hasSurfaceSpecies || hasSurfaceMasters)
     || (hasSurfaceSpecies && hasSurfaceMasters && emptySurfaceSpecies && emptySurfaceMasters)
    )
    {
        surfaceEnabled_ = false;
        surfaceMixturePtr_.clear();
        surfaceMasters_.clear();
        surfaceMasterDensity_.clear();
        return;
    }

    if (!(hasSurfaceSpecies && hasSurfaceMasters))
    {
        FatalErrorInFunction
            << "Surface chemistry requires "
            << "'surfaceSpecies' and 'surfaceMasters'"
            << "in thermoPhysicalProperties."
            << exit(FatalError);
    }

    const wordList surfaceSpeciesNames(thermoDict.subDict("surfaceSpecies").toc());
    surfaceMixturePtr_.reset
    (
        new basicMultiComponentMixture(thermoDict, surfaceSpeciesNames, mesh_)
    );

    surfaceMasters_ = speciesTable(thermoDict.subDict("surfaceMasters").toc());

    bool hasSurfAreaField = mesh_.foundObject<volScalarField>("aSurf");
    aSurf_.assign(static_cast<std::size_t>(mesh_.nCells()),0.0);

    if (hasSurfAreaField)
    {
        const volScalarField& aSurf = mesh_.lookupObject<volScalarField>("aSurf");
        forAll(mesh_.cells(), celli)
        {
            aSurf_[celli] = aSurf[celli];
        }
    }

    surfaceMasterDensity_.setSize(surfaceMasters_.size());

    forAll(surfaceMasters_, i)
    {
        surfaceMasterDensity_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    surfaceMasters_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    nSurf_=0;

    //Get cell index for surfaces

    forAll(mesh_.cells(),celli)
    {
        if (chemistryCellActive(celli) && aSurf_[celli]>1e-3)
        {
            surfaceSolutionIndex_[celli]=nSurf_++;
            surfaceEnabled_=true;
        }
    }

    //equilibrate surface with solution
    forAll(mesh_.cells(),celli)
    {
        if (surfaceMasters_.size()>0)
        {
            // All surface masters currently share the same PHREEQC surface
            // charge and therefore the same reacting-wall area. Use the first
            // master field only to detect cells with boundary surface area.
            const volScalarField::Boundary& Surfbf = surfaceMasterDensity_[0].boundaryField();
            forAll(Surfbf,patchi)
            {
                if (Surfbf[patchi].type() == "reactingWall")
                {
                    const labelList& cellOwner = Surfbf[patchi].patch().faceCells();
                    forAll(Surfbf[patchi],facei)
                    {
                        if
                        (
                            chemistryCellActive(celli)
                         && cellOwner[facei] == celli
                         && surfaceSolutionIndex_[celli] == -1
                        )
                        {
                            surfaceSolutionIndex_[celli]=nSurf_++;
                            surfaceEnabled_=true;
                        }
                    }
                }
            }
        }
    }
}


void phreeqcMixture::initialisePorosityFromField()
{
    //set porosity from eps when available, otherwise use 1.0
    hasPorosityField_ = mesh_.foundObject<volScalarField>("eps");

    if (hasPorosityField_)
    {
        const volScalarField& eps = mesh_.lookupObject<volScalarField>("eps");
        forAll(mesh_.cells(), celli)
        {
            porosity_[celli] = eps[celli];
        }
    }
}


bool phreeqcMixture::chemistryCellActive(const label celli) const
{
    return porosity_[celli] > SMALL
        && porosity_[celli] >= epsChemistryMin_;
}

void phreeqcMixture::printSurfaceConcentration()
{
    if (!surfaceEnabled_ || !surfaceMixturePtr_.valid() || surfConcentration_.empty())
    {
        return;
    }

    const int ncells = mesh_.cells().size();
    forAll(mesh_.cells(), celli)
    {
        forAll(surfaceMixturePtr_().species(), i)
        {
            int surfi = surfaceSolutionIndex_[celli];
            if (surfi > -1)
            {
                Info << surfaceMixturePtr_().species()[i] <<"\n";          
                Info << surfConcentration_[componentSurfaceIndex_[i]*ncells + celli] << "\n";
            }
        }
     }
}


void phreeqcMixture::syncSurfaceStateFromFields()
{
    if (!surfaceEnabled_ || !surfaceMixturePtr_.valid() || surfConcentration_.empty())
    {
        return;
    }

    const int ncells = mesh_.cells().size();
    forAll(surfaceMixturePtr_().species(), i)
    {
        forAll(mesh_.cells(), celli)
        {
            int surfi = surfaceSolutionIndex_[celli];
            if (surfi > -1)
            {
                surfConcentration_
                [
                    componentSurfaceIndex_[i]*ncells + celli
                ] = 0.0;
            }
        }

        volScalarField& Yi = surfaceMixturePtr_().Y(i);
        forAll(mesh_.cells(), celli)
        {
            int surfi = surfaceSolutionIndex_[celli];
            if (surfi > -1 && chemistryCellActive(celli))
            {
                surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] += Yi[celli] / surfArea_[celli]*aSurf_[celli];
            }
        }

        forAll(Yi.boundaryField(), patchi)
        {
            if (Yi.boundaryField()[patchi].type()=="reactingWall")
            {
                const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
                const scalarField& Yfaces = Yi.boundaryField()[patchi];
                const surfaceScalarField& magSf = mesh_.magSf();
                forAll(cellOwner, facei)
                {
                    if (surfArea_[cellOwner[facei]]>0)
                    {
                        surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]] += Yfaces[facei] / surfArea_[cellOwner[facei]] * magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]];
                    }
                }
            }
        }
    }
}

void phreeqcMixture::syncSurfaceStateToFields()
{
    if (!surfaceEnabled_ || !surfaceMixturePtr_.valid())
    {
        return;
    }

    const int ncells = mesh_.cells().size();

    //get surface concentration
    if (nSurfaceSpecies_ > 0)
    {
            forAll(surfaceMixturePtr_().species(), i)
            {
                volScalarField& Yi = surfaceMixturePtr_().Y(i);
                forAll(mesh_.cells(),celli)
                {
                    int surfi = surfaceSolutionIndex_[celli];
                    if (surfi > -1 && chemistryCellActive(celli))
                    {
                        Yi[celli] = surfConcentration_
                        [
                            componentSurfaceIndex_[i]*ncells + celli
                        ]/1000.0;
                    }
                }

                forAll(Yi.boundaryFieldRef(), patchi)
                {
                    if (Yi.boundaryFieldRef()[patchi].type()=="reactingWall")
                    {
                    const labelList& cellOwner = Yi.boundaryFieldRef()[patchi].patch().faceCells();
                    scalarField& Yfaces   = Yi.boundaryFieldRef()[patchi];
                    forAll(cellOwner, facei)
                    {
                        if (saturation_[cellOwner[facei]]>1e-3)
                        {
                            Yfaces[facei] = surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]]/1000;//kmol/m2
                        }
                    }
                }
            }
        }
    }
}

void phreeqcMixture::syncConcentrationFromFields()
{
    const int ncells = mesh_.cells().size();

    forAll(species_, i)
    {
        const volScalarField& Yi = Y_[i];
        const int componentIndex = componentSolutionIndex_[i];

        if (componentIndex < 0)
        {
            continue;
        }

        forAll(mesh_.cells(), celli)
        {
            concentration_[componentIndex*ncells + celli] =
                chemistryCellActive(celli) ? Yi[celli] : 0.0;
        }
    }
}

void phreeqcMixture::printConcentration()
{
    const int ncells = mesh_.cells().size();
    forAll(mesh_.cells(), celli)
    {
        forAll(species_, i)
        {
            const int componentIndex = componentSolutionIndex_[i];

            if (componentIndex < 0)
            {
                continue;
            }

            Info << species_[i] << "\n";
            Info << concentration_[componentIndex*ncells + celli] << "\n";
        }
    }
}

void phreeqcMixture::syncPorosityToModule()
{
    std::vector<double> modulePorosity
    (
        static_cast<std::size_t>(mesh_.cells().size()),
        0.0
    );

    forAll(mesh_.cells(), celli)
    {
        if (chemistryCellActive(celli))
        {
            modulePorosity[celli] = porosity_[celli];
        }
        else
        {
            modulePorosity[celli] = 1.0;
        }
    }

    checkPhreeqcStatus
    (
        "setPorosity",
        phreeqcRM_->setPorosity(modulePorosity.data())
    );
}

// *** initialsie the Phreeqc reactions cells
void phreeqcMixture::initialise()
{
    //get number of cells
    const int ncells = mesh_.cells().size();
    id_=phreeqcRM_->create(ncells, 1);

    Info<< "PhreeqcRM backend: vendored C API, version: 3.3.9-11951"
        << nl;

    gridToChemistry_.assign(static_cast<std::size_t>(ncells), -1);
    chemistryCellCount_ = 0;
    forAll(mesh_.cells(), celli)
    {
        if (chemistryCellActive(celli))
        {
            gridToChemistry_[celli] = chemistryCellCount_++;
        }
    }

    checkPhreeqcStatus
    (
        "createMapping",
        phreeqcRM_->createMapping(gridToChemistry_.data())
    );

    const int backendChemistryCells =
        phreeqcRM_->getChemistryCellCount();
    checkPhreeqcStatus("getChemistryCellCount", backendChemistryCells);

    if (backendChemistryCells != chemistryCellCount_)
    {
        FatalErrorInFunction
            << "PHREEQCRM mapping contains " << backendChemistryCells
            << " cells, but " << chemistryCellCount_ << " were expected."
            << exit(FatalError);
    }

    Info<< "PHREEQCRM fixed chemistry mapping: epsChemistryMin="
        << epsChemistryMin_ << ", activeCells=" << chemistryCellCount_
        << "/" << ncells << nl;

    syncPorosityToModule();

    // Concentrations are represented as mol/L for PHREEQC.
    checkPhreeqcStatus
    (
        "setUnitsSolution",
        phreeqcRM_->setUnitsSolution(2)
    );

    // Save species concentrations for transport coupling.
    checkPhreeqcStatus("setSpeciesSaveOn", phreeqcRM_->setSpeciesSaveOn(true));

    checkPhreeqcStatus
    (
        "loadDatabase",
        phreeqcRM_->loadDatabase("constant/GeoChem.dat")
    );
    checkPhreeqcStatus
    (
        "runFile(constant/phreeqcReactions)",
        phreeqcRM_->runFile(1, 1, 0, "constant/phreeqcReactions")
    );

    if (surfaceEnabled_)
    {
        //solution component list
        std::ostringstream oss;

        //equilibrate surface with solution
        forAll(mesh_.cells(),celli)
        {
            int surfi = surfaceSolutionIndex_[celli];
            if (surfi==-1) continue;

            oss << "SURFACE " << surfi << "\n";
            oss << "-equilibrate with solution " << 0 << "\n";
            forAll(surfaceMasters_, j)
            {
                double area = aSurf_[celli]/1000;
                double mole = surfaceMasterDensity_[j][celli]*aSurf_[celli];
                const volScalarField::Boundary& Surfbf = surfaceMasterDensity_[j].boundaryField();
                forAll(Surfbf,patchi)
                {
                    if (Surfbf[patchi].type() == "reactingWall")
                    {
                        const labelList& cellOwner = Surfbf[patchi].patch().faceCells();
                        const surfaceScalarField& magSf = mesh_.magSf();
                        forAll(Surfbf[patchi],facei)
                        {
                            if (cellOwner[facei]==celli)
                            {
                                mole+=Surfbf[patchi][facei]*magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]];//mol/L
                                area+=magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//m^2/L
                            }
                        }
                    }
                }
                oss << surfaceMasters_[j] << "  " << mole << " " << area  << " 1" << "\n";
            }

            oss << "END\n";
        }

        checkPhreeqcStatus
        (
            "runString(SURFACE templates)",
            phreeqcRM_->runString(0, 1, 0, oss.str().c_str())
        );
    }

    std::vector<int> ic1(static_cast<std::size_t>(7 * ncells), -1);
    forAll(mesh_.cells(), celli)
    {
        if (!chemistryCellActive(celli))
        {
            continue;
        }
        ic1[celli] = 0;               // Solution  i
        ic1[ncells + celli] = -1;           // Equilibrium phases none
        ic1[2 * ncells + celli] = -1;       // Exchange none
        ic1[3 * ncells + celli] = surfaceSolutionIndex_[celli];      // Surface none
        ic1[4 * ncells + celli] = -1;       // Gas phase none
        ic1[5 * ncells + celli] = -1;       // Solid solutions none
        ic1[6 * ncells + celli] = -1;       // Kinetics non
    }

    checkPhreeqcStatus
    (
        "initialPhreeqc2Module",
        phreeqcRM_->initialPhreeqc2Module(ic1.data(), 0, 0)
    );

    checkPhreeqcStatus("findComponents", phreeqcRM_->findComponents());

    //Check what happen when removed
    if (surfaceEnabled_)
    {
        checkPhreeqcStatus
        (
            "runCells(initial surface equilibration)",
            phreeqcRM_->runCells()
        );
    }

    nSolutionSpecies_ = phreeqcRM_->getSpeciesCount();
    checkPhreeqcStatus("getSpeciesCount", nSolutionSpecies_);

    Info<< "nsol:" << nSolutionSpecies_ << "\n";

    std::vector<std::string> solutionComponentNames
    (
        static_cast<std::size_t>(nSolutionSpecies_)
    );

    for (int i = 0; i < nSolutionSpecies_; ++i)
    {
        char nameBuffer[256] = {0};
        phreeqcRM_->getSpeciesName(i, nameBuffer, 256);
        solutionComponentNames[static_cast<std::size_t>(i)] = nameBuffer;
        Info << nameBuffer << endl;
    }

    componentSolutionIndex_.assign(species_.size(), -1);
    forAll(species_, i)
    {
        for (int j = 0; j < nSolutionSpecies_; j++)
        {
            const std::string& component =
                solutionComponentNames[static_cast<std::size_t>(j)];

            if (component == species_[i])
            {
                componentSolutionIndex_[i] = j;
            }
        }
    }

    forAll(species_, i)
    {
        if (componentSolutionIndex_[i] < 0)
        {
            FatalErrorInFunction
                << "Unable to map solution species '" << species_[i]
                << "' to any PhreeqcRM component." << nl
                << ", nsol=" << nSolutionSpecies_
                << exit(FatalError);
        }
    }

    if (nSolutionSpecies_ > 0)
    {
        concentration_.assign
        (
            static_cast<std::size_t>(nSolutionSpecies_ * ncells),
            0.0
        );
        checkPhreeqcStatus
        (
            "getSpeciesConcentrations(initial)",
            phreeqcRM_->getSpeciesConcentrations(concentration_.data())
        );

        //printConcentration();
    }
    else
    {
        concentration_.clear();
    }

    nSurfaceSpecies_ = 0;
    componentSurfaceIndex_.clear();
    surfConcentration_.clear();
    surfArea_.clear();

    // Set water saturation before surface fields can access boundary cells.
    saturation_.assign(static_cast<std::size_t>(ncells), 1.0);
    phreeqcRM_->setSaturation(saturation_.data());

    if (surfaceEnabled_)
    {
        nSurfaceSpecies_ = phreeqcRM_->getSurfaceSpeciesCount();
        checkPhreeqcStatus("getSurfaceSpeciesCount", nSurfaceSpecies_);

        Info<< "nsurf:" << nSurfaceSpecies_ << "\n";

        componentSurfaceIndex_.assign(surfaceMixturePtr_().species().size(), -1);

        if (nSurfaceSpecies_ > 0)
        {
            std::vector<std::string> surfaceComponentNames
            (
                static_cast<std::size_t>(nSurfaceSpecies_)
            );

            for (int i = 0; i < nSurfaceSpecies_; ++i)
            {
                char nameBuffer[256] = {0};
                phreeqcRM_->getSurfaceSpeciesName(i, nameBuffer, 256);
                surfaceComponentNames[static_cast<std::size_t>(i)] = nameBuffer;
            }

            forAll(surfaceMixturePtr_().species(), i)
            {
                for (int j = 0; j < nSurfaceSpecies_; j++)
                {
                    const std::string& component =
                        surfaceComponentNames[static_cast<std::size_t>(j)];

                    if (component == surfaceMixturePtr_().species()[i])
                    {
                        componentSurfaceIndex_[i] = j;
                    }
                }
            }

            forAll(surfaceMixturePtr_().species(), i)
            {
                if (componentSurfaceIndex_[i] < 0)
                {
                    FatalErrorInFunction
                        << "Unable to map surface species '"
                        << surfaceMixturePtr_().species()[i]
                        << "' to any PhreeqcRM surface component." << nl
                        << ", nsurf=" << nSurfaceSpecies_
                        << exit(FatalError);
                }
            }
        }

        if (nSurfaceSpecies_ > 0)
        {
            surfConcentration_.assign
            (
                static_cast<std::size_t>(nSurfaceSpecies_ * ncells),
                0.0
            );
            checkPhreeqcStatus
            (
                "getSurfaceSpeciesConcentrations(initial)",
                phreeqcRM_->getSurfaceSpeciesConcentrations(surfConcentration_.data())
            );
            //printSurfaceConcentration();
            surfArea_.assign(static_cast<std::size_t>(ncells), 0.0);
            phreeqcRM_->getSurfaceArea("Surf", surfArea_.data());
        }
        else
        {
            surfConcentration_.clear();
            surfArea_.clear();
        }
    }

}


phreeqcMixture::phreeqcMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    basicMultiComponentMixture(thermoDict, specieNames, mesh),
    mesh_(mesh),
    surfaceEnabled_(false),
    surfaceMixturePtr_(nullptr),
    surfaceMasters_(),
    surfaceMasterDensity_(),
    hasPorosityField_(false),
    aSurf_(),
    surfaceSolutionIndex_(),
    phreeqcRM_(phreeqcRMAdapter::New()),
    porosity_(static_cast<std::size_t>(mesh.cells().size()), 1.0),
    epsChemistryMin_
    (
        thermoDict.lookupOrDefault<scalar>("epsChemistryMin", 0.0)
    ),
    gridToChemistry_(),
    chemistryCellCount_(0),
    concentration_(),
    componentSolutionIndex_(),
    id_(0),
    nSurf_(0),
    nSolutionSpecies_(0),
    nSurfaceSpecies_(0),
    surfConcentration_(),
    surfArea_(),
    componentSurfaceIndex_(),
    saturation_()
{
    initialisePorosityFromField();
    configureSurfaceMode(thermoDict);
    initialise();
}

phreeqcMixture::phreeqcMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMultiComponentMixture(thermoDict, specieNames, mesh, phaseName),
    mesh_(mesh),
    surfaceEnabled_(false),
    surfaceMixturePtr_(nullptr),
    surfaceMasters_(),
    surfaceMasterDensity_(),
    hasPorosityField_(false),
    aSurf_(),
    surfaceSolutionIndex_(),
    phreeqcRM_(phreeqcRMAdapter::New()),
    porosity_(static_cast<std::size_t>(mesh.cells().size()), 1.0),
    epsChemistryMin_
    (
        thermoDict.lookupOrDefault<scalar>("epsChemistryMin", 0.0)
    ),
    gridToChemistry_(),
    chemistryCellCount_(0),
    concentration_(),
    componentSolutionIndex_(),
    id_(0),
    nSurf_(0),
    nSolutionSpecies_(0),
    nSurfaceSpecies_(0),
    surfConcentration_(),
    surfArea_(),
    componentSurfaceIndex_(),
    saturation_()
{
    initialisePorosityFromField();
    configureSurfaceMode(thermoDict);
    initialise();
}

phreeqcMixture::~phreeqcMixture
(
)
{
    const int destroyStatus = phreeqcRM_->destroy();
    if (destroyStatus != 0)
    {
        WarningInFunction
            << "PhreeqcRM destroy returned non-zero status: "
            << destroyStatus << endl;
    }
}
void phreeqcMixture::correct()
{
    if (concentration_.empty())
    {
        return;
    }

    const int ncells = mesh_.cells().size();

    syncConcentrationFromFields();
    initialisePorosityFromField();
    syncPorosityToModule();
    checkPhreeqcStatus
    (
        "speciesConcentrations2Module(correct)",
        phreeqcRM_->speciesConcentrations2Module(concentration_.data())
    );

    if (surfaceEnabled_ && nSurfaceSpecies_ > 0)
    {
        syncSurfaceStateFromFields();
        checkPhreeqcStatus
        (
            "surfaceSpeciesConcentrations2Module(correct)",
            phreeqcRM_->surfaceSpeciesConcentrations2Module(surfConcentration_.data())
        );
    }


    checkPhreeqcStatus("runCells(correct)", phreeqcRM_->runCells());
    checkPhreeqcStatus
    (
        "getSpeciesConcentrations(correct)",
        phreeqcRM_->getSpeciesConcentrations(concentration_.data())
    );

    if (surfaceEnabled_)
    {
        if (nSurfaceSpecies_ > 0)
        {
            checkPhreeqcStatus
            (
                "getSurfaceSpeciesConcentrations(correct)",
                phreeqcRM_->getSurfaceSpeciesConcentrations(surfConcentration_.data())
            );
        }
    }

    forAll(species_, i)
    {
        volScalarField& Yi = Y_[i];
        const int componentIndex = componentSolutionIndex_[i];

        forAll(mesh_.cells(), celli)
        {
            if (chemistryCellActive(celli) && saturation_[celli]>SMALL)
            {
                Yi[celli] = concentration_[componentIndex * ncells + celli];
            }
        }
    }

    syncSurfaceStateToFields();

}

void phreeqcMixture::setSaturation(const volScalarField& alpha)
{
    forAll(mesh_.cells(), celli)
    {
        if (alpha[celli] >1e-3) saturation_[celli] = alpha[celli];
        else  saturation_[celli]=0;
    }
    phreeqcRM_->setSaturation(saturation_.data());
}

} // End namespace Foam

// ************************************************************************* //
