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

#include "solutionSurfaceMultiComponentMixture.H"
#include "phreeqcRMAdapter.H"
#include "PtrList.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PhreeqcMixture Declaration
\*---------------------------------------------------------------------------*/

class phreeqcMixture
:
    public solutionSurfaceMultiComponentMixture
{
    // Private data
	//initialise phreeqc module 
	void initialise();

	//OpenFOAM mesh
	const fvMesh& mesh_;

	//surface master density fields
	PtrList<volScalarField> surfaceMasterDensity_;

        // chemistry API adapter
        autoPtr<phreeqcRMAdapter> phreeqcRM_;

	//index of phreeqc module
	int id_;

        //number of surfaces in phreeqc
        int nSurf_;

        //component counts
	int nSolutionSpecies_;
	int nSurfaceSpecies_;

	//water saturation
	std::vector<double> saturation_;

	//whether porosity is read from the mesh object registry
	bool hasPorosityField_;

	//porosity passed to PHREEQC
	std::vector<double> porosity_;

        //whether aI is read from the mesh object registry
        bool hasSurfAreaField_;

        //surface area passed to PHREEQC
        std::vector<double> aSurf_;

        //surface solution index
        std::vector<int> surfaceSolutionIndex_;

	//index mapping from solution species to phreeqc component
	std::vector<int> componentSolutionIndex_;
	//index mapping from surface species to phreeqc component
	std::vector<int> componentSurfaceIndex_;
	//concentration
	std::vector<double> concentration_;
	//surface concentration
	std::vector<double> surfConcentration_;
	//surface area
	std::vector<double> surfArea_;


    // Private member functions

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

    // Member functions

	//solve reaction step function
	void correct();

    void setSaturation(const volScalarField& alpha);

        //- Read dictionary
        //void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "phreeqcMixture.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
