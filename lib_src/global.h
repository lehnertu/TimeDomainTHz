/*=========================================================================
 * 
 *  Program:   TimeDomainTHz - THz radiation transport in time-domain
 *
 *  Copyright (c) 2019 U. Lehnert
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#pragma once

/*!
 * \brief defines the natural constants
 * 
 * @author Ulf Lehnert
 * @date 19.11.2019
 * @file global.h
 */

#define ElementaryCharge 1.6021766208e-19       // [As]
#define mecsquared 0.5109989461e6               // electron mass [eV]
#define SpeedOfLight 2.99792458e8               // [m/s]
#define Pi 3.1415926535897932384626433832795029
#define EpsNull 8.854187817e-12                 // vacuum permittivity [As/Vm] 1/c² = e0*µ0
#define MuNull (4*Pi*1e-7)                      // magnetic field constant [Vs/Am]
#define InvRestMass 1.758820025e11              // 1 / m = c² / mc² [m²/s²/eV]

// the electron rest mass in kg can be obtained as
// mecsquared*ElementaryCharge/SpeedOfLight^2

