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

#include "extendedAlgebraicPairGAMGAgglomeration.H"
#include "lduMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedAlgebraicPairGAMGAgglomeration, 0);
    bool extendedAlgebraicPairGAMGAgglomeration::myforward_(true);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        extendedAlgebraicPairGAMGAgglomeration,
        lduMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedAlgebraicPairGAMGAgglomeration::extendedAlgebraicPairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict),
    mymergeLevels_(1)
{
     if (matrix.hasLower())
    {
        myagglomerate
        (
            nCellsInCoarsestLevel_,
            0,
            min(mag(matrix.upper()), mag(matrix.lower())),
            true
        );
    }
    else
    {
        myagglomerate(nCellsInCoarsestLevel_, 0, mag(matrix.upper()), true);
    }      
}


// ************************************************************************* //
