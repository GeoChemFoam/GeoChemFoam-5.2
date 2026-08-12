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
    Foam::multiComponentTransportMixture

Description
    Foam::multiComponentTransportMixture

SourceFiles
    multiComponentTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentTransportMixture_H
#define multiComponentTransportMixture_H

#include "basicMultiComponentTransportMixture.H"
#include "basicMultiComponentMixture.H"
#include "solutionSurfaceMultiComponentMixture.H"
#include "phreeqcMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentTransportMixture Declaration
\*---------------------------------------------------------------------------*/

template <class MixtureType>
class multiComponentTransportMixture
:
    public basicMultiComponentTransportMixture,
    public MixtureType
{
    // Private member functions

        //- Construct as copy (not implemented)
        multiComponentTransportMixture(const multiComponentTransportMixture&);


public:


    // Constructors


        //- Construct from mesh
        multiComponentTransportMixture
        (
            const fvMesh&
        );

    // Member functions

        //- Return the table of species
        const speciesTable& species() const
        {
            return MixtureType::species();
        }

        //- Return the mass-fraction fields
        PtrList<volScalarField>& Y()
        {
            return MixtureType::Y();
        }

        //- Return the const mass-fraction fields
        const PtrList<volScalarField>& Y() const
        {
            return MixtureType::Y();
        }

        //- Return the mass-fraction field for a specie given by index
        volScalarField& Y(const label i)
        {
            return MixtureType::Y(i);
        }

        //- Return the const mass-fraction field for a specie given by index
        const volScalarField& Y(const label i) const
        {
            return MixtureType::Y(i);
        }

        //- Update properties
        void correct()
        {
            MixtureType::correct();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "multiComponentTransportMixture.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
