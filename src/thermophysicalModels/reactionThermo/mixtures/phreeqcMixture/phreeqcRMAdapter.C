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

#include "phreeqcRMAdapter.H"
#include "RM_interface_C.h"
#include "error.H"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace Foam
{

phreeqcRMAdapter::phreeqcRMAdapter()
:
    id_(-1),
    nxyz_(0)
{}

phreeqcRMAdapter::~phreeqcRMAdapter()
{}

    int phreeqcRMAdapter::create(int nxyz, int nthreads)
    {
        nxyz_ = nxyz;
        id_ = RM_Create(nxyz, nthreads);
        return id_;
    }

    int phreeqcRMAdapter::destroy()
    {
        if (id_ < 0)
        {
            return 0;
        }

        const int status = RM_Destroy(id_);
        id_ = -1;
        nxyz_ = 0;
        return status;
    }

    int phreeqcRMAdapter::createMapping(int* gridToChemistry)
    {
        return RM_CreateMapping(id_, gridToChemistry);
    }

    int phreeqcRMAdapter::getChemistryCellCount()
    {
        return RM_GetChemistryCellCount(id_);
    }

    int phreeqcRMAdapter::setPorosity(double* porosity)
    {
        return RM_SetPorosity(id_, porosity);
    }

    int phreeqcRMAdapter::setUnitsSolution(int units)
    {
        return RM_SetUnitsSolution(id_, units);
    }

    int phreeqcRMAdapter::setSpeciesSaveOn(bool saveOn)
    {
        return RM_SetSpeciesSaveOn(id_, saveOn ? 1 : 0);
    }

    int phreeqcRMAdapter::loadDatabase(const char* database)
    {
        return RM_LoadDatabase(id_, database);
    }

    int phreeqcRMAdapter::runFile(int workers, int initialPhreeqc, int utility, const char* chemistryName)
    {
        return RM_RunFile(id_, workers, initialPhreeqc, utility, chemistryName);
    }

    int phreeqcRMAdapter::runString(int workers, int initialPhreeqc, int utility, const char* inputString)
    {
        return RM_RunString(id_, workers, initialPhreeqc, utility, inputString);
    }

    int phreeqcRMAdapter::initialPhreeqc2Module(int* ic1, int* ic2, double* fraction1)
    {
        return RM_InitialPhreeqc2Module(id_, ic1, ic2, fraction1);
    }

    int phreeqcRMAdapter::runCells()
    {
        return RM_RunCells(id_);
    }

    int phreeqcRMAdapter::findComponents()
    {
        return RM_FindComponents(id_);
    }

    int phreeqcRMAdapter::getSpeciesCount()
    {
        return RM_GetSpeciesCount(id_);
    }

    int phreeqcRMAdapter::getSpeciesName(int index, char* name, int length)
    {
        return RM_GetSpeciesName(id_, index, name, length);
    }

    int phreeqcRMAdapter::getSurfaceSpeciesCount()
    {
        return RM_GetSurfaceSpeciesCount(id_);
    }

    int phreeqcRMAdapter::getSurfaceSpeciesName(int index, char* name, int length)
    {
        return RM_GetSurfaceSpeciesName(id_, index, name, length);
    }

    int phreeqcRMAdapter::getSpeciesConcentrations(double* concentrations)
    {
        return RM_GetSpeciesConcentrations(id_, concentrations);
    }

    int phreeqcRMAdapter::getSurfaceSpeciesConcentrations(double* concentrations)
    {
        return RM_GetSurfaceSpeciesConcentrations(id_, concentrations);
    }

    int phreeqcRMAdapter::getSurfaceArea(const char* surfaceName, double* area)
    {
        return RM_GetSurfaceArea(id_, surfaceName, area);
    }

    int phreeqcRMAdapter::setSaturation(double* saturation)
    {
        return RM_SetSaturation(id_, saturation);
    }

    int phreeqcRMAdapter::speciesConcentrations2Module(double* concentrations)
    {
        return RM_SpeciesConcentrations2Module(id_, concentrations);
    }

    int phreeqcRMAdapter::surfaceSpeciesConcentrations2Module(double* concentrations)
    {
        return RM_SurfaceSpeciesConcentrations2Module(id_, concentrations);
    }

    int phreeqcRMAdapter::decodeError(int errorCode) const
    {
        return RM_DecodeError(id_, errorCode);
    }

    std::string phreeqcRMAdapter::getErrorString() const
    {
        const int len = RM_GetErrorStringLength(id_);
        if (len <= 0)
        {
            return std::string();
        }

        std::string error
        (
            static_cast<std::size_t>(len) + 1u,
            '\0'
        );

        if (RM_GetErrorString(id_, &error[0], len + 1) != 0)
        {
            return std::string();
        }

        error.resize(std::strlen(error.c_str()));
        return error;
    }

autoPtr<phreeqcRMAdapter> phreeqcRMAdapter::New()
{
    return autoPtr<phreeqcRMAdapter>(new phreeqcRMAdapter());
}

} // End namespace Foam

// ************************************************************************* //
