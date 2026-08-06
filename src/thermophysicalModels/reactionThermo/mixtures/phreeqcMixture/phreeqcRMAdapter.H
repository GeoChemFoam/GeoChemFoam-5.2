/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef phreeqcRMAdapter_H
#define phreeqcRMAdapter_H

#include "autoPtr.H"
#include "word.H"
#include <string>

namespace Foam
{

class phreeqcRMAdapter
{
protected:

    int id_;
    int nxyz_;

public:

    phreeqcRMAdapter();

    virtual ~phreeqcRMAdapter();

    int create(int nxyz, int nthreads);

    int destroy();

    int setPorosity(double* porosity);

    int setUnitsSolution(int units);

    int setSpeciesSaveOn(bool saveOn);

    int loadDatabase(const char* database);

    int runFile(int workers, int initialPhreeqc, int utility, const char* chemistryName);

    int runString(int workers, int initialPhreeqc, int utility, const char* inputString);

    int initialPhreeqc2Module(int* ic1, int* ic2, double* fraction1);

    int runCells();

    int findComponents();

    int getSpeciesCount();

    int getSpeciesName(int index, char* name, int length);

    int getSurfaceSpeciesCount();

    int getSurfaceSpeciesName(int index, char* name, int length);

    int getSpeciesConcentrations(double* concentrations);

    int getSurfaceSpeciesConcentrations(double* concentrations);

    int getSurfaceArea(const char* surfaceName, double* area);

    int setSaturation(double* saturation);

    int speciesConcentrations2Module(double* concentrations);

    int surfaceSpeciesConcentrations2Module(double* concentrations);

    static autoPtr<phreeqcRMAdapter> New();
};

} // End namespace Foam

#endif

// ************************************************************************* //
