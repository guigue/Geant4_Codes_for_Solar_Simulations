//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
// Adapted from:
// 
//  S.Larsson and J. Generowicz.
//
//    *************************************
//    *                                   *
//    *    PurgMagTabulatedField3D.cc     *
//    *                                   *
//    *************************************
//
//

#include "Solar_TabulatedField.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

namespace{
  G4Mutex mySolarTabulatedFieldLock = G4MUTEX_INITIALIZER;
}

Solar_TabulatedField::Solar_TabulatedField(const char* filename)
  :invertX(false),invertY(false),invertZ(false)
{

  double lenUnit= km;
  double fieldUnit= tesla;

  G4cout << "\n-----------------------------------------------------------"
         << "\n      Magnetic field"
         << "\n-----------------------------------------------------------";

  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << endl;

  //
  //This is a thread-local class and we have to avoid that all workers open the
  //file at the same time
  G4AutoLock lock(&mySolarTabulatedFieldLock);

  ifstream file( filename ); // Open the file for reading.

  if (!file.is_open())
    {
      G4ExceptionDescription ed;
      ed << "Could not open input file " << filename << G4endl;
      G4Exception("SolarTabulatedField::SolarTabulatedField",
      "pugmag001",FatalException,ed);
    }

  // Ignore first blank line
  char buffer[256];
  file.getline(buffer,256);

  // Read table dimensions
  file >> nx >> ny >> nz; // Note dodgy order

  G4cout << "  [ Number of values x,y,z: "
   << nx << " " << ny << " " << nz << " ] "
   << endl;

  // Set up storage space for table
  xField.resize( nx );
  yField.resize( nx );
  zField.resize( nx );
  int ix, iy, iz;
  for (ix=0; ix<nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }

  // Ignore other header information
  // The first line whose second character is '0' is considered to
  // be the last line of the header.
  do {
    file.getline(buffer,256);
  } while ( buffer[1]!='0');
  // Change y and x order of reading to read Python dipole April19
  // Read in the data
  double xval,yval,zval,bx,by,bz;
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
        file >> xval >> yval >> zval >> bx >> by >> bz;
        if ( ix==0 && iy==0 && iz==0 ) {
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx * fieldUnit;
        yField[ix][iy][iz] = by * fieldUnit;
        zField[ix][iy][iz] = bz * fieldUnit;
      }
    }
  }

  file.close();

  lock.unlock();

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  G4cout << "\n ---> ... done reading " << endl;

  // G4cout << " Read values of field from file " << filename << endl;
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
   << "\n ---> Min values x,y,z: "
   << minx/km << " " << miny/km << " " << minz/km << " km "
   << "\n ---> Max values x,y,z: "
   << maxx/km << " " << maxy/km << " " << maxz/km << " km " << endl;
  // << "\n ---> The field will be offset by " << zOffset/cm << " cm " << endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;}
  if (maxy < miny) {swap(maxy,miny); invertY = true;}
  if (maxz < minz) {swap(maxz,minz); invertZ = true;}
  G4cout << "\nAfter reordering if neccesary"
   << "\n ---> Min values x,y,z: "
   << minx/km << " " << miny/km << " " << minz/km << " km "
   << " \n ---> Max values x,y,z: "
   << maxx/km << " " << maxy/km << " " << maxz/km << " km ";

  dx = maxx - minx;
  epsilon_x = dx / (nx-1) ;
  dy = maxy - miny;
  epsilon_y = dy / (ny-1) ;
  dz = maxz - minz;
  epsilon_z = dz / (nz-1) ;

  G4cout << "\n ---> Dif values x,y,z (range): "
   << dx/km << " " << dy/km << " " << dz/km << " km in z "
   << "\n-----------------------------------------------------------" << endl;

}

void Solar_TabulatedField::GetFieldValue(const double point[4],
              double *Bfield ) const
{

  double x = point[0];
  double y = point[1];
  double z = point[2]; //+ fZoffset;

  // Check that the point is within the defined region
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy &&
       z>=minz && z<=maxz )
    {

    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy;
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;

    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));

    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);

    // Full 3-dimensional version
    Bfield[0] =
      xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

    Bfield[1] =
      yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

    Bfield[2] =
      zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //dBx/dx
    Bfield[3] = (Bfield[0]  -
		 (xField[xindex][yindex  ][zindex  ]  *  (1-ylocal) * (1-zlocal) +
		  xField[xindex][yindex  ][zindex+1]  *  (1-ylocal) *    zlocal  +
		  xField[xindex][yindex+1][zindex  ]  *   ylocal    * (1-zlocal) +
		  xField[xindex][yindex+1][zindex+1]  *   ylocal    *    zlocal
		  )) / (xlocal * epsilon_x);

    // dBx/dy
    Bfield[4] = (Bfield[0] -
		 (xField[xindex  ][yindex  ][zindex  ]  * (1-xlocal) * (1-zlocal) +
		  xField[xindex  ][yindex  ][zindex+1]  * (1-xlocal) *    zlocal  +
		  xField[xindex+1][yindex  ][zindex  ]  *    xlocal  * (1-zlocal) +
		  xField[xindex+1][yindex  ][zindex+1]  *    xlocal  *    zlocal
		  )) / (ylocal * epsilon_y);

    //dBx/dz
    Bfield[5] = (Bfield[0] -
		 (xField[xindex  ][yindex  ][zindex  ]  * (1-xlocal) * (1-ylocal)  +
		  xField[xindex  ][yindex+1][zindex  ]  * (1-xlocal) *    ylocal   +
		  xField[xindex+1][yindex  ][zindex  ]  *    xlocal  * (1-ylocal)  +
		  xField[xindex+1][yindex+1][zindex  ]  *    xlocal  *    ylocal
		  )) / (zlocal * epsilon_z);

    // dBy/dx
    Bfield[6] = (Bfield[1] -
		 (yField[xindex  ][yindex  ][zindex  ]  * (1-ylocal) * (1-zlocal) +
		  yField[xindex  ][yindex  ][zindex+1] * (1-ylocal) *    zlocal  +
		  yField[xindex  ][yindex+1][zindex  ]  *    ylocal  * (1-zlocal) +
		  yField[xindex  ][yindex+1][zindex+1] *    ylocal  *    zlocal
		  )) / (xlocal * epsilon_x);
    // dBy/dy
    Bfield[7] = (Bfield[1] -
		 (yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-zlocal) +
		  yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) *    zlocal  +
		  yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-zlocal) +
		  yField[xindex+1][yindex  ][zindex+1] *    xlocal  *    zlocal
		  )) / (ylocal * epsilon_y);

    // dBy/dz
    Bfield[8] = (Bfield[1] -
		 (yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) +
		  yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  +
		  yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) +
		  yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal
		  )) / (zlocal * epsilon_z);

    // dBz/dx
    Bfield[9] = (Bfield[2] -
		 (zField[xindex  ][yindex  ][zindex  ] * (1-ylocal) * (1-zlocal) +
		  zField[xindex  ][yindex  ][zindex+1] * (1-ylocal) *    zlocal  +
		  zField[xindex  ][yindex+1][zindex  ] *    ylocal  * (1-zlocal) +
		  zField[xindex  ][yindex+1][zindex+1] *    ylocal  *    zlocal
		  )) / (xlocal * epsilon_x);

    // dBz/dy
    Bfield[10] = (Bfield[2] -
		  (zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-zlocal) +
		   zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) *    zlocal  +
		   zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-zlocal) +
		   zField[xindex+1][yindex  ][zindex+1] *    xlocal  *    zlocal
		   )) / (ylocal * epsilon_y);

    // dBz/dz
    Bfield[11] = (Bfield[2] -
		  (zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) +
		   zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  +
		   zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) +
		   zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal
		   )) / (zlocal * epsilon_z);

    } else {

    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
    Bfield[3] = 0.0;
    Bfield[4] = 0.0;
    Bfield[5] = 0.0;
    Bfield[6] = 0.0;
    Bfield[7] = 0.0;
    Bfield[8] = 0.0;
    Bfield[9] = 0.0;
    Bfield[10] = 0.0;
    Bfield[11] = 0.0;
  }

  //G4cout << ' ' << endl;

}
