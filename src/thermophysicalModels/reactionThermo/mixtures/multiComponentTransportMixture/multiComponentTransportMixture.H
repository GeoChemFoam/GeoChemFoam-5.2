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

#include "IOdictionary.H"
#include "basicMultiComponentMixture.H"
#include "solutionSurfaceMultiComponentMixture.H"
#include "phreeqcMixture.H"
#include "reactiveMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentTransportMixture Declaration
\*---------------------------------------------------------------------------*/

template <class MixtureType>
class multiComponentTransportMixture
:
    public IOdictionary,
    public MixtureType
{
    // Private data

        //- Species data
        PtrList<dimensionedScalar> DY_;

        //- Flag saying whether species i has a tensorial DTensor field
        boolList hasTensorDY_;

        //- Species tensorial diffusion/dispersion tensor
        //  DTensorY_[i] is only allocated if hasTensorDY_[i] == true
        PtrList<volTensorField> DTensorY_;

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

        //- Read dictionary
        //void read(const dictionary&);
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
