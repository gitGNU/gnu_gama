/*  
    GNU Gama --- Geodesy and Mapping C++ library 
    Copyright (C) 1999, 2003  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
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
 *  $Id: version.cpp,v 1.4 2003/05/11 12:32:25 cepek Exp $
 */


#include <gnu_gama/version.h>

namespace GNU_gama {

const char* GNU_gama_version  = "1.7.04-pre";

const char* GNU_gama_compiler =
              #if   defined (__GNUC__)
              "GNU g++"             // g++ 3.0.4 and 2.95.4
              #elif defined (__BORLANDC__) && (__linux__)
              "kylix-bc++"          // 5.7
              #elif defined (__BORLANDC__)
              "win32-borland"       // 5.6
              #elif defined (_MSC_VER)
              "win32-msvc"          // 1300
              #else
              #error GNU_gama - has not been tested with your compiler
              #endif
              ;
}


/* GNU Gama uses James Clark's parser Expat (v.1.1) for XML data processing
 *
 *    Expat is subject to the Mozilla Public License Version 1.1. 
 *    Alternatively you may use expat under the GNU General Public
 *    License instead.
 *
 *              ftp://ftp.jclark.com/pub/xml/expat.zip
 *
 * Scripts for compiling GNU Gama and linking the program gama-local expect
 * Expat library to be in the same directory as GNU Gama

=============================================================================

1.7.04 2003-05-11

   - classes BaseParser and Dataparser moved from namespace GaMaLib to
     namespace GNU_gama.


1.7.03 2003-03-29

   - first draft of implementation of (g3) classes for computing
     derivatives of Observations by Parameters

   - files gamalib/version.[h|cpp] (ChangeLog) moved to directory gnu_gama

   - bash scripts for generating makefiles replaced by C++ program
     scripts/gnu_gama_dep.cpp

     * for building makefiles under win32 see scripts/win32-makefiles.bat  

     * gamalib.a renamed to libgama.a

   - adjustments needed by MSVC compiler (in files gnu_gama/list.h and
     gamalib/cluster.h)


1.7.02 2003-03-13

   - added template class PointBase (gnu_gama/pointbase.h)

   - added template class template <class T> class List<T*> (file
     gnu_gama/list.h); shall replace all pointer lists used in Gama

   - removed bug in template BlockDiagonal (sparse/sbdiagonal.h)
     repoted by Jan Bilek

   - removed bug in SparserMatrix (gnu_gama/sparse/smatrix.h) reported
     by Jan Bilek; not enough memory was allocated for rptr array in
     SparseMatrix(const SparseMatrix* sm)

   - added missing 'typenames' into gnu_gama/obsdata.h


1.7.01 2003-02-28 

   - observation data structures (formerly defined in gamalib/cluster.*
     and gamalib/local/gamadata.*) rewritten as template functions 
     and moved into gnu_gama/obsdata.h

   - directory `sparse' moved from `gamalib' to `gnu_gama', added
     template class SparseVec (added directory <gnu_gama/...> for
     processing in scripts/gamalaib_dep.cpp)


1.7.00 2003-02-16

     #####################################################################
     #                                                                   #
     #  unstable version 1.7 (main trunk at CVS savannah.gnu.org)        #
     #                                                                   #
     #####################################################################

   - for the new development branch (g3) is introduced new namespace
     GNU_gama; this name is going to replace old name GaMaLib
     (similarly directory gnu_gama is going to replace directory
     gamalib)

   - added output encoding option into gama-local program. Output encoding
     is handled by class GNU_gama::OutStream

   - removed bugs in functions `underline' and `set_width' (bad text
     length in utf-8 encoding)


1.6.00 2003-02-01

     #####################################################################
     #                                                                   #
     #  stable version 1.6 (branch gama-1-6 at CVS savannah.gnu.org)     #
     #                                                                   #
     #####################################################################

   - at CVS savannah.gnu.org stable version 1.6 (adjustment in local
     coordinate system, program gama-local) has been started as a
     separate development branch with label gama-1-6; development of
     unstable version 1.7 (aimed to adjustment in geocentric
     coordinates, program gama-g3) is going to continue in the CVS
     main trunk.

    - Jan Pytel: reduction of observations with nonzero heights of
      instrument and/or target in gama-local


1.5.08 2002-01-18

    - Jan Pytel: bug in zangle.h  (missing test for coordinates xy)

    - added data object "AdjInputDataObject" into DataParser

    - first draft of classes DataParser, AdjInputData and Adj


1.5.07 2002-12-18

    - added missing declarations "using namespace std;" in various files
      (apart from g++ 2.95.4 tested with bcc32 5.6 and cl 13.00.9466).

         Note: Makefiles for win platform are not converted using unix2dos

    - call to std:fmod() in g2d_helper.h replaced with expresion "x%2"


1.5.06 2002-12-14

    - added reductions for basic observation types (Jan Pytel)

    - Dutch translation of output texts produced by gama-local (John Dedrum)


1.5.05 2002-11-21

    - removed bug in text output of fixed coordinates in gama-local
      reported by Carl Verheyen (in the case of "inconsistent
      coordinate system" internally changed sign of y coordinates was
      not transformed back in the output).  File: fixed_points.h

    - gmatvec-0.9.21 (removed forgotten calls to fabs() )

    - some minor changes in scripts/GaMaLib_archive to reflect new
      archive file names in format 'gama-xx.yy.zz.tar.gz'


1.5.04 2002-10-25

    - updated documantation (added references to 'gama-local')

    - removed command line program 'gama' 


1.5.03 2002-10-24

    - tar archive file names changed from 'gamalib-xx.yy.zz.tar.gz' to
      'gama-xx.yy.zz.tar.gz'

    - class Point renamed to LocalPoint (file local/point.h renamed to
      local/lpoint.h)

    - command line program 'gama' renamed to 'gama-local' (all previous
      scripts needed for building 'gama' are left in 1.5.03 archive
      tar and will be removed in a next version)

    - general parts of GKFparser moved into the new class BaseParser.
      Class BaseParser is public base class to GKFparser and a new
      GaMa XML parser DataParser

    - added a test if dos2unix utility is available into
      scripts/Build_GaMa (Jan Pytel)


1.5.02 2002-09-29

   - support for the Finish language

   - added scripts/kylix.sh for compilation GaMa with bc++ compilere
     from Kylix3 Open Edition under Linux


1.5.01 2002-09-15

   - new version of BandMat<> class (gmatvec 0.9.20)

   - initial version of new class BlockDiagonal (in gamalib/sparse)
     for symmetric block diagonal matrix (planned for new adjustment
     model)

   - expat source included into gamalib archive as a convenience for
     GaMa users

   - initial draft of the Ellipsoid class and a helper program
     ellipsoids_xml.cpp for generating files ellipsoids.[h|cpp] from
     xml/ellipsoids.xml

   - merged changes from 1.4.01:

     -- changes concerning heights of instrument  and/or reflector
     -- removed unused public data member `tag' from the class Cluster
        (in files cluster.h and observation.[h|cpp])

   - removed bug reported by Jan Bilek in utf8_cp1250(char *buf); file
     gamalib/xml/encoding.cpp:

         310c310
         <       iso_8859_2_unicode((int*)tab);
         ---
         >       cp1250_unicode((int*)tab);


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











