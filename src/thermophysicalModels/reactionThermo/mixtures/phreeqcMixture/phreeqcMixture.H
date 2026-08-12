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

Class
    Foam::phreeqcMixture

Description
    Foam::phreeqcMixture

SourceFiles
    phreeqcMixture.C

\*---------------------------------------------------------------------------*/

#ifndef phreeqcMixture_H
#define phreeqcMixture_H

#include "basicMultiComponentMixture.H"
#include "phreeqcRMAdapter.H"
#include "PtrList.H"
#include "speciesTable.H"
#include <vector>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PhreeqcMixture Declaration
\*---------------------------------------------------------------------------*/

class phreeqcMixture
:
    public basicMultiComponentMixture
{
    // Private data

        //- OpenFOAM mesh
        const fvMesh& mesh_;

        //- Optional cell-based surface mixture state
        bool surfaceEnabled_;

        //surface species mixture
        autoPtr<basicMultiComponentMixture> surfaceMixturePtr_;
        speciesTable surfaceMasters_;
        PtrList<volScalarField> surfaceMasterDensity_;

        //whether porosity is read from the mesh object registry
        bool hasPorosityField_;

        //surface area - 0 if on present
        std::vector<double> aSurf_;

        //surface solution index
        std::vector<int> surfaceSolutionIndex_;

        //- Chemistry API adapter
        autoPtr<phreeqcRMAdapter> phreeqcRM_;

        //- Physical mobile-water porosity used by PHREEQC
        std::vector<double> porosity_;

        //- Minimum fixed eps for inclusion in the PHREEQCRM mapping
        scalar epsChemistryMin_;

        //- OpenFOAM-cell to PHREEQCRM chemistry-cell mapping
        std::vector<int> gridToChemistry_;

        //- Number of cells included in the fixed PHREEQCRM mapping
        label chemistryCellCount_;

        //- Solution species concentrations, ordered by PHREEQC component index
        std::vector<double> concentration_;

        //- Index mapping from OpenFOAM solution species to PHREEQC component
        std::vector<int> componentSolutionIndex_;

	//index of phreeqc module
	int id_;

        //number of surfaces in phreeqc
        int nSurf_;

        //- Cached component count from active backend
        int nSolutionSpecies_;
        int nSurfaceSpecies_;

        //- Surface-species concentrations, ordered by PHREEQC surface index
        std::vector<double> surfConcentration_;

        //- Surface-species area per liter water for the PHREEQC backend
        std::vector<double> surfArea_;

        //- Index mapping from OpenFOAM surface species to PHREEQC component
        std::vector<int> componentSurfaceIndex_;

	//water saturation
	std::vector<double> saturation_;

    // Private Member Functions

        //- Fail fast on PHREEQC backend errors instead of continuing into UB
        void checkPhreeqcStatus(const char* operation, int status) const;

        //- Configure optional cell-based surface support from thermoDict
        void configureSurfaceMode(const dictionary& thermoDict);

        //- Populate the porosity cache from field if present
        void initialisePorosityFromField();

        //- print Phreeqc concentration
        void printConcentration();

        //- Copy the current OpenFOAM solution fields into the PHREEQC buffer
        void syncConcentrationFromFields();

        //- print Phreeqc surface concentration
        void printSurfaceConcentration();

        //- Copy the current OpenFOAM surface fields into the PHREEQC buffer
        void syncSurfaceStateFromFields();

        //- Push the PHREEQC surface state back into OpenFOAM surface fields
        void syncSurfaceStateToFields();

        //- Initialise the PHREEQC module for aqueous-only DBS chemistry
        void initialise();

        //- Push a backend-safe porosity vector into PHREEQC
        void syncPorosityToModule();

        //- Return true when a cell belongs to the fixed chemistry mapping
        bool chemistryCellActive(const label celli) const;

        //- Populate the current cellwise surface area scale from aSurf
        void syncSurfaceAreaFromFields();

        //- Construct as copy (not implemented)
        phreeqcMixture(const phreeqcMixture&);


public:

    // Constructors

        //- Construct from mesh and object registry
        phreeqcMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&
        );

        //- Construct from mesh and object registry and phase name
        phreeqcMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh&,
            const word&
        );

    // Destructor
    virtual ~phreeqcMixture();

    // Member Functions

        //- Solve reaction step
        void correct();

    void setSaturation(const volScalarField& alpha);


};

} // End namespace Foam

#endif

// ************************************************************************* //
