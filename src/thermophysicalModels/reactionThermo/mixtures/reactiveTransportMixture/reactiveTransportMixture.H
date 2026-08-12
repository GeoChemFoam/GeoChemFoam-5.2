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
    Foam::reactiveTransportMixture

Description
    Reactive transport mixture wrapper around a multiComponentTransportMixture.

SourceFiles
    reactiveTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactiveTransportMixture_H
#define reactiveTransportMixture_H

#include "basicReactiveTransportMixture.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class reactiveTransportMixture Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class reactiveTransportMixture
:
    public basicReactiveTransportMixture,
    public MixtureType
{

protected:

        //- Index mapping to index of 1st SI component
        std::vector<int> component1Index_;

        //- Index mapping to index of 2nd SI component
        std::vector<int> component2Index_;


private:

        //- Construct as copy (not implemented)
        reactiveTransportMixture(const reactiveTransportMixture&);


public:

    //- Runtime type information
    TypeName("reactiveTransportMixture");


    // Constructors

        //- Construct from mesh
        reactiveTransportMixture(const fvMesh&);


    //- Destructor
    virtual ~reactiveTransportMixture();


    // Member functions

        //- Correct mixture and kinetic reaction state
        void correct();

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

        void setSaturation(const volScalarField& alpha)
        {
            // Nothing to do for the kinetic wrapper.
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
