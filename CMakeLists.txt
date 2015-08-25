cmake_minimum_required(VERSION 2.8)
project(GNU_gama)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=\\\"1.16X\\\" ")

include_directories(lib)

set(GAMA_SOURCE_FILES
##    bin/gama-g3.cpp
##    bin/gama-local-xml2sql.cpp
##    bin/gama-local-xml2txt.cpp
##    bin/gama-local.cpp
#    lib/expat/gennmtab/gennmtab.c    #-DGNU_gama_expat_1_1
#    lib/expat/sample/elements.c
##    lib/expat/xmlparse/hashtable.c
##    lib/expat/xmlparse/hashtable.h
##    lib/expat/xmlparse/xmlparse.c
##    lib/expat/xmlparse/xmlparse.h
##    lib/expat/xmltok/asciitab.h
#    lib/expat/xmltok/dllmain.c
##    lib/expat/xmltok/iasciitab.h
##    lib/expat/xmltok/latin1tab.h
##    lib/expat/xmltok/nametab.h
##    lib/expat/xmltok/utf8tab.h
##    lib/expat/xmltok/xmldef.h
##    lib/expat/xmltok/xmlrole.c
##    lib/expat/xmltok/xmlrole.h
##    lib/expat/xmltok/xmltok.c
##    lib/expat/xmltok/xmltok.h
#    lib/expat/xmltok/xmltok_impl.c
#    lib/expat/xmltok/xmltok_impl.h
#    lib/expat/xmltok/xmltok_ns.c
##    lib/expat/xmlwf/codepage.c
##    lib/expat/xmlwf/codepage.h
#    lib/expat/xmlwf/filemap.h
#    lib/expat/xmlwf/readfilemap.c
#    lib/expat/xmlwf/unixfilemap.c
#    lib/expat/xmlwf/win32filemap.c
#    lib/expat/xmlwf/xmlfile.c
#    lib/expat/xmlwf/xmlfile.h
#    lib/expat/xmlwf/xmltchar.h
#    lib/expat/xmlwf/xmlwf.c
    lib/gnu_gama/adj/adj.cpp
    lib/gnu_gama/adj/adj.h
    lib/gnu_gama/adj/adj_base.h
    lib/gnu_gama/adj/adj_basefull.h
    lib/gnu_gama/adj/adj_basesparse.h
    lib/gnu_gama/adj/adj_chol.h
    lib/gnu_gama/adj/adj_chol_implementation.h
    lib/gnu_gama/adj/adj_envelope.h
    lib/gnu_gama/adj/adj_envelope_implementation.h
    lib/gnu_gama/adj/adj_gso.h
    lib/gnu_gama/adj/adj_svd.h
    lib/gnu_gama/adj/envelope.h
    lib/gnu_gama/adj/homogenization.h
    lib/gnu_gama/g3/g3_adjres.cpp
    lib/gnu_gama/g3/g3_adjres.h
    lib/gnu_gama/g3/g3_cluster.cpp
    lib/gnu_gama/g3/g3_cluster.h
    lib/gnu_gama/g3/g3_model.cpp
    lib/gnu_gama/g3/g3_model.h
    lib/gnu_gama/g3/g3_model_init.cpp
    lib/gnu_gama/g3/g3_model_linearization.cpp
    lib/gnu_gama/g3/g3_model_revision.cpp
    lib/gnu_gama/g3/g3_model_write_xml_adjustment_results.cpp
    lib/gnu_gama/g3/g3_observation.h
    lib/gnu_gama/g3/g3_parameter.h
    lib/gnu_gama/g3/g3_point.cpp
    lib/gnu_gama/g3/g3_point.h
    lib/gnu_gama/g3/g3_write_observation_xml.cpp
    lib/gnu_gama/g3/g3_write_observation_xml.h
    lib/gnu_gama/local/acord/approx_heights.cpp
    lib/gnu_gama/local/acord/approx_heights.h
    lib/gnu_gama/local/acord/approx_vectors.cpp
    lib/gnu_gama/local/acord/approx_vectors.h
    lib/gnu_gama/local/acord/reduce_observations.cpp
    lib/gnu_gama/local/acord/reduce_observations.h
    lib/gnu_gama/local/acord/reduce_to_ellipsoid.cpp
    lib/gnu_gama/local/acord/reduce_to_ellipsoid.h
    lib/gnu_gama/local/median/g2d_cogo.cpp
    lib/gnu_gama/local/median/g2d_cogo.h
    lib/gnu_gama/local/median/g2d_coordinates.cpp
    lib/gnu_gama/local/median/g2d_coordinates.h
    lib/gnu_gama/local/median/g2d_exception.h
    lib/gnu_gama/local/median/g2d_helper.cpp
    lib/gnu_gama/local/median/g2d_helper.h
    lib/gnu_gama/local/median/g2d_point.cpp
    lib/gnu_gama/local/median/g2d_point.h
    lib/gnu_gama/local/pobs/format.cpp
    lib/gnu_gama/local/pobs/format.h
    lib/gnu_gama/local/pobs/write.cpp
    lib/gnu_gama/local/results/text/adjusted_observations.h
    lib/gnu_gama/local/results/text/adjusted_unknowns.h
    lib/gnu_gama/local/results/text/approximate_coordinates.h
    lib/gnu_gama/local/results/text/error_ellipses.h
    lib/gnu_gama/local/results/text/fixed_points.h
    lib/gnu_gama/local/results/text/general_parameters.h
    lib/gnu_gama/local/results/text/network_description.h
    lib/gnu_gama/local/results/text/outlying_abs_terms.h
    lib/gnu_gama/local/results/text/reduced_observations.h
    lib/gnu_gama/local/results/text/reduced_observations_to_ellipsoid.h
    lib/gnu_gama/local/results/text/residuals_observations.h
    lib/gnu_gama/local/results/text/underline.cpp
    lib/gnu_gama/local/results/text/underline.h
    lib/gnu_gama/local/acord.cpp
    lib/gnu_gama/local/acord.h
    lib/gnu_gama/local/angobs.h
    lib/gnu_gama/local/bearing.cpp
    lib/gnu_gama/local/bearing.h
    lib/gnu_gama/local/bpoint.h
    lib/gnu_gama/local/c_api.h
    lib/gnu_gama/local/cluster.h
    lib/gnu_gama/local/display_observation_visitor.cpp
    lib/gnu_gama/local/display_observation_visitor.h
    lib/gnu_gama/local/exception.h
    lib/gnu_gama/local/float.h
    lib/gnu_gama/local/gamadata.cpp
    lib/gnu_gama/local/gamadata.h
    lib/gnu_gama/local/html.cpp
    lib/gnu_gama/local/html.h
    lib/gnu_gama/local/itstream.h
    lib/gnu_gama/local/language.cpp
    lib/gnu_gama/local/language.h
    lib/gnu_gama/local/lcoords.h
    lib/gnu_gama/local/linearization.h
    lib/gnu_gama/local/local_linearization_visitor.cpp
    lib/gnu_gama/local/local_linearization_visitor.h
    lib/gnu_gama/local/local_revision_visitor.cpp
    lib/gnu_gama/local/local_revision_visitor.h
    lib/gnu_gama/local/localnetwork2sql.cpp
    lib/gnu_gama/local/localnetwork2sql.h
    lib/gnu_gama/local/lpoint.h
    lib/gnu_gama/local/matvec.h
    lib/gnu_gama/local/medianf.h
    lib/gnu_gama/local/network.cpp
    lib/gnu_gama/local/network.h
    lib/gnu_gama/local/observation.cpp
    lib/gnu_gama/local/observation.h
    lib/gnu_gama/local/orientation.cpp
    lib/gnu_gama/local/orientation.h
    lib/gnu_gama/local/pointid.cpp
    lib/gnu_gama/local/pointid.h
    lib/gnu_gama/local/readsabw.h
    lib/gnu_gama/local/revision.h
    lib/gnu_gama/local/skipcomm.cpp
    lib/gnu_gama/local/skipcomm.h
    lib/gnu_gama/local/sqlitereader.cpp
    lib/gnu_gama/local/sqlitereader.h
    lib/gnu_gama/local/svg.cpp
    lib/gnu_gama/local/svg.h
    lib/gnu_gama/local/test_linearization_visitor.cpp
    lib/gnu_gama/local/test_linearization_visitor.h
    lib/gnu_gama/local/writevisitor.h
    lib/gnu_gama/local/xmlerror.cpp
    lib/gnu_gama/local/xmlerror.h
    lib/gnu_gama/sparse/intlist.h
    lib/gnu_gama/sparse/sbdiagonal.h
    lib/gnu_gama/sparse/smatrix.h
    lib/gnu_gama/sparse/smatrix_graph.h
    lib/gnu_gama/sparse/smatrix_graph_connected.h
    lib/gnu_gama/sparse/smatrix_ordering.h
    lib/gnu_gama/sparse/svector.h
    lib/gnu_gama/xml/baseparser.cpp
    lib/gnu_gama/xml/baseparser.h
    lib/gnu_gama/xml/dataobject.cpp
    lib/gnu_gama/xml/dataobject.h
    lib/gnu_gama/xml/dataparser.cpp
    lib/gnu_gama/xml/dataparser.h
    lib/gnu_gama/xml/dataparser_adj.cpp
    lib/gnu_gama/xml/dataparser_g3.cpp
    lib/gnu_gama/xml/dataparser_g3adj.cpp
    lib/gnu_gama/xml/encoding.cpp
    lib/gnu_gama/xml/encoding.h
    lib/gnu_gama/xml/encoding_cp1251.cpp
    lib/gnu_gama/xml/encoding_unknown_handler.cpp
    lib/gnu_gama/xml/gkfparser.cpp
    lib/gnu_gama/xml/gkfparser.h
    lib/gnu_gama/xml/htmlparser.cpp
    lib/gnu_gama/xml/htmlparser.h
    lib/gnu_gama/xml/localnetwork_adjustment_results.cpp
    lib/gnu_gama/xml/localnetwork_adjustment_results.h
    lib/gnu_gama/xml/localnetwork_adjustment_results_data.cpp
    lib/gnu_gama/xml/localnetwork_adjustment_results_data.h
    lib/gnu_gama/xml/localnetworksql.cpp
    lib/gnu_gama/xml/localnetworksql.h
    lib/gnu_gama/xml/localnetworkxml.cpp
    lib/gnu_gama/xml/localnetworkxml.h
    lib/gnu_gama/e3.cpp
    lib/gnu_gama/e3.h
    lib/gnu_gama/ellipsoid.cpp
    lib/gnu_gama/ellipsoid.h
    lib/gnu_gama/ellipsoids.cpp
    lib/gnu_gama/ellipsoids.h
    lib/gnu_gama/exception.h
    lib/gnu_gama/gon2deg.cpp
    lib/gnu_gama/gon2deg.h
    lib/gnu_gama/intfloat.h
    lib/gnu_gama/latlong.cpp
    lib/gnu_gama/latlong.h
    lib/gnu_gama/list.h
    lib/gnu_gama/model.h
    lib/gnu_gama/movetofront.h
    lib/gnu_gama/obsdata.h
    lib/gnu_gama/outstream.cpp
    lib/gnu_gama/outstream.h
    lib/gnu_gama/pointbase.h
    lib/gnu_gama/radian.h
    lib/gnu_gama/rand.cpp
    lib/gnu_gama/rand.h
    lib/gnu_gama/statan.cpp
    lib/gnu_gama/statan.h
    lib/gnu_gama/utf8.cpp
    lib/gnu_gama/utf8.h
    lib/gnu_gama/version.cpp
    lib/gnu_gama/version.h
    lib/gnu_gama/visitor.h
    lib/gnu_gama/xml_expat.h
    lib/gnu_gama/xsd.h
    lib/matvec/array.h
    lib/matvec/bandmat.h
    lib/matvec/choldec.h
    lib/matvec/covmat.h
    lib/matvec/gso.h
    lib/matvec/hilbert.h
    lib/matvec/inderr.h
    lib/matvec/jacobian.h
    lib/matvec/mat.h
    lib/matvec/matbase.h
    lib/matvec/matvec.h
    lib/matvec/matvecbase.h
    lib/matvec/memrep.h
    lib/matvec/pinv.h
    lib/matvec/sortvec.h
    lib/matvec/svd.h
    lib/matvec/symmat.h
    lib/matvec/transmat.h
    lib/matvec/transvec.h
    lib/matvec/unsigned.h
    lib/matvec/vec.h
    lib/matvec/vecbase.h
##    scripts/ellipsoids_xml.cpp
##    scripts/gnu_gama_dep.cpp
##    scripts/slovnikar.cpp
##    tests/gama-local/scripts/check_algorithms.cpp
##    tests/gama-local/scripts/check_equivalents.cpp
##    tests/gama-local/scripts/check_html.cpp
##    tests/gama-local/scripts/check_xml_coordinates.cpp
##    tests/gama-local/scripts/check_xml_results.cpp
##    tests/gama-local/scripts/check_xml_xml.cpp
##    tests/gama-local/scripts/check_xyz.cpp
##    tests/gama-local/scripts/check_xyz.h
##    tests/gama-local/scripts/compare_xml_adjustment.cpp
##    tests/gama-local/scripts/compare_xml_adjustment.h
##    tests/gama-local/scripts/sqlite_init_db.cpp
)

add_library(gama STATIC ${GAMA_SOURCE_FILES})

include(FindEXPAT)
find_package(EXPAT REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

add_executable(gama-local bin/gama-local.cpp)
target_link_libraries(gama-local gama ${EXPAT_LIBRARIES})

#add_executable(gama-local-xml2txt bin/gama-local-xml2txt.cpp)
#target_link_libraries(gama-local gama ${EXPAT_LIBRARIES})
