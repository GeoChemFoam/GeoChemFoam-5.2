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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::phreeqcMixture::initialise()
{
    //get number of cells
    int ncells = mesh_.cells().size();
    id_=phreeqcRM_->create(ncells, 1);

    //set porosity from eps when available, otherwise use 1.0
    hasPorosityField_ = mesh_.foundObject<volScalarField>("eps");
    porosity_.assign(static_cast<std::size_t>(ncells), 1.0);

    if (hasPorosityField_)
    {
        const volScalarField& eps = mesh_.lookupObject<volScalarField>("eps");
        forAll(mesh_.cells(), celli)
        {
            porosity_[celli] = eps[celli];
        }
    }

    hasSurfAreaField_ = mesh_.foundObject<volScalarField>("aSurf");
    aSurf_.assign(static_cast<std::size_t>(ncells),0.0);

    if (hasSurfAreaField_)
    {
        const volScalarField& aSurf = mesh_.lookupObject<volScalarField>("aSurf");
        forAll(mesh_.cells(), celli)
        {
            aSurf_[celli] = aSurf[celli];
        }
    }

    phreeqcRM_->setPorosity(porosity_.data());
    
    //concentration=mol/L
    phreeqcRM_->setUnitsSolution(2);

    //Save species for multi-species transport
    phreeqcRM_->setSpeciesSaveOn(true);

    //load phreeqc database
    phreeqcRM_->loadDatabase("constant/GeoChem.dat");

    //load reaction
    phreeqcRM_->runFile(1, 1, 0, "constant/phreeqcReactions");

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


    //solution component list
    std::ostringstream oss;


    //init surface and solution
    forAll(mesh_.cells(),celli)
    {
        oss << "COPY solution 0 " << celli  << "\n";
        oss << "END" << "\n";
    }

    surfaceSolutionIndex_.assign(static_cast<std::size_t>(ncells), -1);

    nSurf_=0;
    //Get cell index for surfaces
    forAll(mesh_.cells(),celli)
    {
        if (aSurf_[celli]>1e-3)
        {
            surfaceSolutionIndex_[celli]=nSurf_++;
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
                        if (cellOwner[facei]==celli && surfaceSolutionIndex_[celli]==-1)
                        {
                            surfaceSolutionIndex_[celli]=nSurf_++;
                        }
                    }
                }
            }
        }
    }

    //equilibrate surface with solution
    forAll(mesh_.cells(),celli)
    {
        int surfi = surfaceSolutionIndex_[celli];
        if (surfi==-1) continue;

        oss << "SURFACE " << surfi << "\n";
        oss << "-equilibrate with solution " << celli << "\n";
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
        oss << "END" << "\n";
    }

    Info << oss.str() << endl;
    //run phreeqc keywords
    phreeqcRM_->runString(0, 1, 0, oss.str().c_str());

    //init Phreeqc worker module
    std::vector<int> ic1(static_cast<std::size_t>(7 * ncells), -1);
    forAll(mesh_.cells(),celli)
    {
        ic1[celli] = celli;               // Solution  i
        ic1[ncells + celli] = -1;      // Equilibrium phases none
        ic1[2 * ncells + celli] = -1;       // Exchange none
        ic1[3 * ncells + celli] = surfaceSolutionIndex_[celli];      // Surface none
        ic1[4 * ncells + celli] = -1;      // Gas phase none
        ic1[5 * ncells + celli] = -1;      // Solid solutions none
        ic1[6 * ncells + celli] = -1;      // Kinetics none
    }

    phreeqcRM_->initialPhreeqc2Module(ic1.data(), 0, 0);

    //find components
    phreeqcRM_->findComponents();
    
    //get number of solution species
    nSolutionSpecies_ = phreeqcRM_->getSpeciesCount();

    //display number of silution species
    Info << "nsol:" << nSolutionSpecies_ << "\n";

    std::vector<std::string> solutionComponentNames
    (
        static_cast<std::size_t>(nSolutionSpecies_)
    );
    for (int i = 0; i < nSolutionSpecies_; ++i)
    {
        char nameBuffer[256] = {0};
        phreeqcRM_->getSpeciesName(i, nameBuffer, 256);
        solutionComponentNames[static_cast<std::size_t>(i)] = nameBuffer;
    }

    //get number of surface species
    nSurfaceSpecies_ = phreeqcRM_->getSurfaceSpeciesCount();

    //display number of surface species
    Info << "nsurf:" << nSurfaceSpecies_ << "\n";

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
    
    //save component index map
    componentSolutionIndex_.assign(species_.size(), -1);
    componentSurfaceIndex_.assign(surfaceSpecies_.size(), -1);
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

    forAll(surfaceSpecies_, i)
    {
        for (int j = 0; j < nSurfaceSpecies_; j++)
        {
            const std::string& component =
                surfaceComponentNames[static_cast<std::size_t>(j)];
            if (component == surfaceSpecies_[i])
            {
                componentSurfaceIndex_[i] = j;
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

    //concentration, adsorption and surface potential
    if (nSolutionSpecies_>0)
    {
        concentration_.assign
        (
            static_cast<std::size_t>(nSolutionSpecies_ * ncells),
            0.0
        );
    }
    else
    {
        concentration_.clear();
    }
    if (nSurfaceSpecies_>0) 
    {
        surfConcentration_.assign
        (
            static_cast<std::size_t>(nSurfaceSpecies_ * ncells),
            0.0
        );
        surfArea_.assign(static_cast<std::size_t>(ncells), 0.0);
    }
    else
    {
        surfConcentration_.clear();
        surfArea_.clear();
    }
    
    //Run phreeqc
    //phreeqcRM_->runCells();

    if (nSolutionSpecies_ > 0)
    {
        phreeqcRM_->getSpeciesConcentrations(concentration_.data());
    }
    if (nSurfaceSpecies_ > 0)
    {
        phreeqcRM_->getSurfaceSpeciesConcentrations(surfConcentration_.data());
        phreeqcRM_->getSurfaceArea("Surf", surfArea_.data());
    }

    //set water saturation for Phreeqc module
    saturation_.assign(static_cast<std::size_t>(ncells), 1.0);
    phreeqcRM_->setSaturation(saturation_.data());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phreeqcMixture::phreeqcMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    solutionSurfaceMultiComponentMixture(thermoDict, specieNames, mesh),
    mesh_(mesh),
    surfaceMasterDensity_(),
    phreeqcRM_(phreeqcRMAdapter::New()),
    nSurf_(0),
    nSolutionSpecies_(0),
    nSurfaceSpecies_(0),
    saturation_(),
    hasPorosityField_(false),
    porosity_(),
    hasSurfAreaField_(false),
    aSurf_(),
    surfaceSolutionIndex_(),
    componentSolutionIndex_(),
    componentSurfaceIndex_(),
    concentration_(),
    surfConcentration_(),
    surfArea_()
{
    initialise();
}

Foam::phreeqcMixture::phreeqcMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    solutionSurfaceMultiComponentMixture(thermoDict, specieNames, mesh, phaseName),
    mesh_(mesh),
    surfaceMasterDensity_(),
    phreeqcRM_(phreeqcRMAdapter::New()),
    nSurf_(0),
    nSolutionSpecies_(0),
    nSurfaceSpecies_(0),
    saturation_(),
    hasPorosityField_(false),
    porosity_(),
    hasSurfAreaField_(false),
    aSurf_(),
    surfaceSolutionIndex_(),
    componentSolutionIndex_(),
    componentSurfaceIndex_(),
    concentration_(),
    surfConcentration_(),
    surfArea_()
{
    initialise();
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::phreeqcMixture::~phreeqcMixture
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::phreeqcMixture::correct()
{

    //get number of cells
    int ncells = mesh_.cells().size();

    if (hasPorosityField_)
    {
        const volScalarField& eps = mesh_.lookupObject<volScalarField>("eps");
        forAll(mesh_.cells(), celli)
        {
            porosity_[celli] = eps[celli];
        }
        phreeqcRM_->setPorosity(porosity_.data());
    }

    //get concentration after transport
    forAll(species_, i)
    {
        volScalarField& Yi = Y_[i];
        forAll(mesh_.cells(), celli)
        {
            concentration_[componentSolutionIndex_[i] * ncells + celli] = Yi[celli];//mol/L
        }
    }


    //get surface concentration after transport
    if (nSurfaceSpecies_ > 0)
    {
            forAll(surfaceSpecies_, i)
            {
                forAll(mesh_.cells(), celli)
                {
                    int surfi = surfaceSolutionIndex_[celli];
                    if (surfi>-1) surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] = 0.0;
                }

                volScalarField& Yi = surfaceMixture_.Y(i);
                forAll(mesh_.cells(),celli)
                {
                    int surfi = surfaceSolutionIndex_[celli];
                    if (surfi>-1) surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] += Yi[celli] / surfArea_[celli]*aSurf_[celli];
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

    //set concentration for Phreeqc module
    phreeqcRM_->speciesConcentrations2Module(concentration_.data());
    if (nSurfaceSpecies_ > 0)
    {
        phreeqcRM_->surfaceSpeciesConcentrations2Module(surfConcentration_.data());
    }

    //Run phreeqc
    phreeqcRM_->runCells();

    //get concentration after reactions
    phreeqcRM_->getSpeciesConcentrations(concentration_.data());
    if (nSurfaceSpecies_ > 0)
    {
        phreeqcRM_->getSurfaceSpeciesConcentrations(surfConcentration_.data());
    }

    //get reaction rate and new vector composition
    forAll(species_, i)
    {
        volScalarField& Yi = Y_[i];
        forAll(mesh_.cells(), celli)
        {
            if (saturation_[celli]>1e-3)
            {
                Yi[celli]  = concentration_[componentSolutionIndex_[i] * ncells + celli];//mol/m3
            }
        }
    }

    //get surface concentration
    if (nSurfaceSpecies_ > 0)
    {
            forAll(surfaceSpecies_, i)
            {
                volScalarField& Yi = surfaceMixture_.Y(i);
                forAll(mesh_.cells(),celli)
                {
                    int surfi = surfaceSolutionIndex_[celli];
                    if (surfi>-1) Yi[celli]=surfConcentration_[componentSurfaceIndex_[i] * ncells + celli]/1000;
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

void Foam::phreeqcMixture::setSaturation(const volScalarField& alpha)
{
    forAll(mesh_.cells(), celli)
    {
        if (alpha[celli] >1e-3) saturation_[celli] = alpha[celli];
        else  saturation_[celli]=0;
    }
    phreeqcRM_->setSaturation(saturation_.data());
}
// ************************************************************************* //
