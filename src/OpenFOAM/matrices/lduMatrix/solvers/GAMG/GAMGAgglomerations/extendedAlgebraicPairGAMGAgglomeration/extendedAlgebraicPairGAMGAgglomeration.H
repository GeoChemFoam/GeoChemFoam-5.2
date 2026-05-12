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
    Foam::extendedAlgebraicPairGAMGAgglomeration

Description
    Agglomerate using the extended pair algorithm.

SourceFiles
    extendedAlgebraicPairGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef extendedAlgebraicPairGAMGAgglomeration_H
#define extendedAlgebraicPairGAMGAgglomeration_H

#include "pairGAMGAgglomeration.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class extendedAlgebraicPairGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class extendedAlgebraicPairGAMGAgglomeration
:
    public pairGAMGAgglomeration
{

    // Private data

        //- Number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
        label mymergeLevels_;

        //- Direction of cell loop for the current level
        static bool myforward_;
        
public:

    //- Runtime type information
    TypeName("extendedAlgebraicPair");


    // Constructors

        //- Construct given mesh and controls
        extendedAlgebraicPairGAMGAgglomeration
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );
        
        //- Calculate and return agglomeration
        static tmp<labelField> myagglomerate
        (
            label& nCoarseCells,
            const lduAddressing& fineMatrixAddressing,
            const scalarField& faceWeights
        );

        //- Agglomerate from a starting level. Starting level is usually 0
        //- (initial mesh) but sometimes >0 (restarting after processor
        //- agglomeration)
        virtual void myagglomerate
        (
            const label nCellsInCoarsestLevel,
            const label startLevel,
            const scalarField& startFaceWeights,
            const bool doProcessorAgglomerate = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
