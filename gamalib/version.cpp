/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 1999  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU GaMa / GaMaLib C++ Library.
    
    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 *  $Id: version.cpp,v 1.26 2002/06/11 19:19:37 cepek Exp $
 */


#include <gamalib/version.h>

namespace GaMaLib {

const char* GaMaLib_version  = "1.5.01-pre";

const char* GaMaLib_compiler =
              #ifdef __GNUC__
              "GNU g++"             // g++ 2.95.2
              #elif defined __BORLANDC__
              "win32-borland"       // 5.5
              #elif defined _MSC_VER
              "win32-msvc"          // 6.0
              #else
              #error GaMaLib - has not been tested with your compiler
              #endif
              ;
}


/* GaMaLib uses James Clark's parser Expat (v.1.1) for XML data processing
 *
 *    Expat is subject to the Mozilla Public License Version 1.1. 
 *    Alternatively you may use expat under the GNU General Public
 *    License instead.
 *
 *              ftp://ftp.jclark.com/pub/xml/expat.zip
 *
 * Scripts for compiling GaMaLib and linking the program GaMa expect
 * Expat library to be in the same directory as GaMaLib

=============================================================================

1.5.01-pre

     - initial draft of the Ellipsoid class and a helper program
       ellipsoids_xml.cpp for generating files ellipsoids.[h|cpp] from
       xml/ellipsoids.xml


1.5.00 2002-06-11

     #####################################################################
     #                                                                   #
     #  unstable version 1.5 (main trunk at CVS savannah.gnu.org)        #
     #                                                                   #
     #####################################################################

   - unstable version 1.5 is aimed at adjustment of geodetic network
     in geocentric coordinate system. Apart from this comment, version
     1.5.00 is identical to 1.4.00.


1.4.00 2002-06-11

     #####################################################################
     #                                                                   #
     #  stable version 1.4 (branch gama-1-4 at CVS savannah.gnu.org)     #
     #                                                                   #
     #####################################################################

   - at CVS savannah.gnu.org stable version 1.4 (adjustment in local
     coordinate system) has been started as a separate branch with label
     gama-1-4; development of unstable version 1.5 (aimed to adjustment
     in geocentric coordinates) is going to continue in the CVS main
     trunk.


1.3.39 2002-06-09

   - various bugs reported by g++ 3.0.4 (Jan Pytel); files
     observation.h, g2d_coordinates.cpp, statan.cpp, write.cpp.
     
   - a bug in helper program gamalib_dep.cpp (Jan Pytel)

   - some minor changes in makefiles; makefiles and language.[h|cpp] 
     added to CVS repository

   - directory gamalib/gamaxml renamed to gamalib/xml


1.3.38 2002-05-29

   - class Point : added functions x_0(), y_0(), z_0() and set_xyz_0()
     used for correct output of dx/dy/dz corrections in the case of
     iterated adjustment (output of coordinates corrections for the
     initial approximation and not to the last iteration).

   - removed directory gamalib/doc and files scripts/build-doc-html and
     build scripts/gamalib2xhtml*

     automaticaly generated html files for browsing gamalib sources were
     neither maintained nor used


1.3.37 2002-05-24

   - new attributes bs="..." and fs="..." in the definition of <angle />

     Chuck Ghilani pointed out that XML description of angles was not
     natural; thus attribute `to' was renamed to `bs' (backsight) and
     `rs' was renamed to `fs' (foresight)

     Former attributes `to' and `rs' of the tag <angle /> was left as
     an undocumented feature in the class  GKFparser

   - function Angle::rs() was renamed to Angle::fs() and added new
     function Angle::bs() which is a synonym for Angle::to()


1.3.36 2002-05-15

   - in the class GKFparser implicitly all covariance matrices are
     checked for positive-definiteness. This test can be suppressed
     if the function GKFparser::check_covariances(false) is called.
     
     This dark feature was introduced to GKFparser by Jan Pytel. If
     you use it you go on your own risk. You have been warned!


1.3.35 2002-04-06

   - acord.cpp inproved by Jan Pytel (removed unnecessary iterations)


1.3.34 2002-04-05

   - approximate coordinates are computed in iterations if necessary
     (gamalib/local/acord.cpp); max. number of iterations is set to 5.


1.3.33 2002-04-03

   - new scripts for generating Makefiles (inspired by qmake Makefiles
     and Qt project files). All platform dependent environment variables
     are defined in a single file scripts/platforms.defs.


1.3.32 2002-03-19

   - gmatvec 0.9.17 (a bug in BandMat::cholDec(); Jan Pytel: Cholesky
     decomposition failed if the first element of the matrix was zero)
   
   - gmatvec 0.9.16 (added iterators to Mat<> / Vec<> classes)


1.3.31 2001-12-21

   - in scripts/make-linux2borland.sed added parameters needed by
     Rocinante and a minor change in scripts/Build_GaMa (Jan Pytel)

   - gmatvec 0.9.15: a bug in SymMat::invert() reported by Leos Mervart;
     in the case of dimension==1 inversion was computed twice
     and thus nothing happened (missing return statement)

   - a bug in computation of approximate coordinates reported by Jiri
     Vesely; if there were two neighbouring directions with identical
     targets, gama reported exception and stopped (the exception was
     thrown by the Angle constructor).

     From 1.3.31 we allow for angles left and right targets to be
     identical, surely this is not an realistic case but it's useful
     in g2d_point.h:67 and should not cause a trouble elsewhere


1.3.30  2001-12-07

   - a bug in ObservationData::deepCopy found by Jan Pytel (covariance
     matrices were not copied)

       diff -r1.1 gamadata.cpp
       59a60
       >       current->covariance_matrix = (*ci)->covariance_matrix;


1.3.29  2001-12-02

   - GNU GaMa sources installed at Savannah CVS repository

   - old change logs 

               from    version 1.3.28  2001-11-30
               back to         0.1.00  1998-01-19 

     are recorded in the file ChangeLog.1

*/











