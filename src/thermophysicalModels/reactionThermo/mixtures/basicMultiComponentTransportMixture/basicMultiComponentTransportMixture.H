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
    Foam::basicMultiComponentTransportMixture

Description
    Base transport class for multicomponent mixtures. This class owns
    transport data but does not own mixture composition fields.

SourceFiles
    basicMultiComponentTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicMultiComponentTransportMixture_H
#define basicMultiComponentTransportMixture_H

#include "IOdictionary.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class basicMultiComponentTransportMixture Declaration
\*---------------------------------------------------------------------------*/

class basicMultiComponentTransportMixture
:
    public IOdictionary
{

protected:

    // Protected data

        //- Species data
        PtrList<dimensionedScalar> DY_;

        //- Flag saying whether species i has a tensorial DTensor field
        boolList hasTensorDY_;

        //- Species tensorial diffusion/dispersion tensor
        //  DTensorY_[i] is only allocated if hasTensorDY_[i] == true
        PtrList<volTensorField> DTensorY_;


private:

        //- Construct as copy (not implemented)
        basicMultiComponentTransportMixture
        (
            const basicMultiComponentTransportMixture&
        );


public:

    // Constructors

        //- Construct from mesh
        basicMultiComponentTransportMixture(const fvMesh&);


    //- Destructor
    virtual ~basicMultiComponentTransportMixture();


    // Member functions

        //- Return the table of species
        virtual const speciesTable& species() const = 0;

        //- Return the mass-fraction fields
        virtual PtrList<volScalarField>& Y() = 0;

        //- Return the const mass-fraction fields
        virtual const PtrList<volScalarField>& Y() const = 0;

        //- Return the mass-fraction field for a specie given by index
        virtual volScalarField& Y(const label i) = 0;

        //- Return the const mass-fraction field for a specie given by index
        virtual const volScalarField& Y(const label i) const = 0;

        //- Return the species molecular diffusion list
        const PtrList<dimensionedScalar>& DY() const
        {
            return DY_;
        }

        //- Return the species molecular diffusion
        dimensionedScalar DY(const label i) const
        {
            return DY_[i];
        }

        //- Return the species tensorial D list
        const PtrList<volTensorField>& DTensorY() const
        {
            return DTensorY_;
        }

        //- Return the tensorial DTensor for species i
        const volTensorField& DTensorY(const label i) const
        {
            return DTensorY_[i];
        }

        //- Return whether species i has a tensorial DTensor
        bool hasTensorDY(const label i) const
        {
            return hasTensorDY_[i];
        }

        //- Return the full list of tensorial DTensor flags
        const boolList& hasTensorDY() const
        {
            return hasTensorDY_;
        }

        //- Update properties
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
