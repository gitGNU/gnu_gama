#include <gamalib/c_api.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* $Id: capi_demo.c,v 1.1 2001/12/07 12:33:41 cepek Exp $ */

/* This simple program demonstrates usage of C API functions to
 * GaMa/GaMaLib C++ classes and functions. See gama-api.texinfo for
 * more detailed description */

int main(int argc, char *argv[])
{
  const char* exception;
  char*       network_description;
  char        buf[1024];

  if (argc < 2 || argc > 3)
    {
      printf("\nusage: %s  input_file [output_file]\n\n", argv[0]);
      return 1;
    }


  /* "do { } while(0); if() { }" simulates C++ try/catch blocks  */
  do 
    {
      void* local_network, *file, *gkf_parser;

      /* initalization of C API (language: 0 == English) */
      Cgama_init(0);

      /* create new LocalNetwork object (algorithm SVD) */
      local_network = Cgama_LocalNetwork_svd();
      if (Cgama_exception()) break;

      /* setting output file name */
      if (argc == 3)
        {
          strcpy(buf, argv[2]);
        }
      else
        {
          strcpy(buf, argv[1]);
          strcat(buf, "-");
          strcat(buf, Cgama_gamalib_version());
          strcat(buf, "-");
          strcat(buf, Cgama_LocalNetwork_algorithm(local_network));
        }
      
      /* C API output file constructor */
      file = Cgama_output_file(local_network, buf);
      
      Cgama_of_string(file, "\n----  GNU GaMa C API Demo  --------------  ");
      Cgama_of_string(file, argv[1]);
      Cgama_of_string(file, " ----\n\n");

      /* read and parse XML input data */
      do {
        /* GKF parser constructor */
        gkf_parser = Cgama_GKF_parser(local_network);
        if (Cgama_exception()) break;

        {
          size_t  n;
          FILE   *inp;

          inp = fopen(argv[1], "r");
          while((n = fread(buf, 1, 1024, inp)))
            {
              Cgama_GKF_parser_parse(gkf_parser, buf, n, 0);
            }
          Cgama_GKF_parser_parse(gkf_parser, buf, 0, 1);
        }

        /* parameters for statistical analysis */
        Cgama_LocalNetwork_set_apriori_m0
          (local_network, Cgama_GKF_parser_apriori_m0(gkf_parser) );
        Cgama_LocalNetwork_set_conf_pr
          (local_network, Cgama_GKF_parser_conf_pr(gkf_parser) );
        Cgama_LocalNetwork_set_tol_abs
          (local_network, Cgama_GKF_parser_tol_abs(gkf_parser) );
        Cgama_LocalNetwork_set_type_refsd
          (local_network, Cgama_GKF_parser_m0_apriori(gkf_parser) );
        
        /* network description as a C string */
        network_description = Cgama_GKF_parser_description(gkf_parser);

        /* GKF parser destructor */
        Cgama_GKF_parser_dtor(gkf_parser);

      }              
      while (0);
      if (Cgama_exception()) break;

      /* check if we have something to do */
      if (Cgama_LocalNetwork_PointData_empty(local_network))
        {
          Cgama_of_string(file, "No points available\n");
          Cgama_output_file_close(file);
          return 1;
        }
      if (Cgama_LocalNetwork_ObservationData_empty(local_network))
        {
          Cgama_of_string(file, "No observations available\n");
          Cgama_output_file_close(file);
          return 1;
        }

      /* compute and print approximate coordinates */
      Cgama_of_approximate_coordinates(file);
      if (Cgama_exception()) break;

      /* check if we have something to do */
      if (Cgama_LocalNetwork_sum_points(local_network)   == 0 ||
          Cgama_LocalNetwork_sum_unknowns(local_network) == 0 )
        {
          Cgama_of_string(file, "No network points defined\n");
          Cgama_output_file_close(file);
          return 1;
        }

      /* outlying absolute terms */
      if (Cgama_LocalNetwork_huge_abs_terms(local_network))
        {
          Cgama_of_outlying_abs_terms(file);
          Cgama_LocalNetwork_remove_huge_abs_terms(local_network);
          Cgama_of_string(file, 
                          "Observations with huge abs. terms removed\n\n");
        }

      /* print user's description of the network*/
      if (network_description)
        {
          Cgama_of_network_description(file, network_description);
          free(network_description);
        }

      /* now check if the given network can be adjusted and print all
         adjustment results */
      if (Cgama_of_general_parameters(file))
        {
          int iteration = 0;
          do
            {
              if (++iteration > 1)
                {
                  Cgama_of_string(file,
                                  "Approximate coordinates corrected\n");
                  Cgama_of_string(file,
                                  "*********************************\n\n");
                  Cgama_LocalNetwork_refine_approx(local_network);
                  Cgama_of_general_parameters(file);
                }
              Cgama_of_fixed_points(file);
              Cgama_of_adjusted_unknowns(file);
            }
          while (Cgama_of_test_linearization(file) && iteration < 3);
        
        }
      
      Cgama_of_error_ellipses(file);
      Cgama_of_adjusted_observations(file);
      Cgama_of_residuals_observations(file);


      /* output file destructor */
      Cgama_output_file_close(file);

      /* delete LocalNetwork */
      Cgama_LocalNetwork_dtor(local_network);
      if (Cgama_exception()) break;

    }
  while(0);
  if ((exception = Cgama_exception()))
    {
      if (Cgama_exception_unknown()) 
        exception = "unknown exception";
      
      printf("program ended with error: %s \n", exception);
      {
        int n = Cgama_exception_GKF_parser_line();
        if (n)
          {
            printf("XML input file %s :  line %d\n", argv[1], n);
          }
      }
    }
  
  return 0;
}










