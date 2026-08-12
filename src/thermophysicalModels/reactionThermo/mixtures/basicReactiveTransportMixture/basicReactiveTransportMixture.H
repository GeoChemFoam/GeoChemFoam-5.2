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
    Foam::basicReactiveTransportMixture

Description
    Base reactive-transport state and logic independent of the underlying
    mixture type.

SourceFiles
    basicReactiveTransportMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicReactiveTransportMixture_H
#define basicReactiveTransportMixture_H

#include "autoPtr.H"
#include "basicMultiComponentTransportMixture.H"
#include "dimensionedScalar.H"
#include "PtrList.H"
#include "runTimeSelectionTables.H"
#include "speciesTable.H"
#include "volFields.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class basicReactiveTransportMixture Declaration
\*---------------------------------------------------------------------------*/

class basicReactiveTransportMixture
:
    virtual public basicMultiComponentTransportMixture
{

protected:

    // Protected data

        wordList kineticPhase_;
        wordList kineticPhaseReaction_;

        PtrList<dimensionedScalar> Mws_;
        PtrList<dimensionedScalar> rhos_;
        PtrList<dimensionedScalar> Ke_;
        PtrList<dimensionedScalar> k0_;
        PtrList<scalarList> scoeff_;
        PtrList<labelList> kiSpeciesIndex_;
        PtrList<PtrList<dimensionedScalar> > ki_;

        PtrList<volScalarField> Ws_;
        PtrList<volScalarField> Omega_;
        PtrList<volScalarField> R_;

    // Protected member functions


private:

        //- Construct as copy (not implemented)
        basicReactiveTransportMixture(const basicReactiveTransportMixture&);


public:

    //- Runtime type information
    TypeName("basicReactiveTransportMixture");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicReactiveTransportMixture,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        basicReactiveTransportMixture(const fvMesh&);


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<basicReactiveTransportMixture> New(const fvMesh& mesh);


    //- Destructor
    virtual ~basicReactiveTransportMixture();


    // Member functions

        const wordList& kineticPhases()
        {
            return kineticPhase_;
        }

        const wordList& kineticPhaseReactions()
        {
            return kineticPhaseReaction_;
        }

        const dimensionedScalar Mws(label i)
        {
            return Mws_[i];
        }

        const dimensionedScalar rhos(label i)
        {
            return rhos_[i];
        }

        const volScalarField& Ws(label i)
        {
            return Ws_[i];
        }

        const volScalarField& omega(label i)
        {
            return Omega_[i];
        }

        const scalarList& scoeff(label i)
        {
            return scoeff_[i];
        }

        const volScalarField& Rate(label i)
        {
            return R_[i];
        }

        //- Update properties
        virtual void correct()=0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
