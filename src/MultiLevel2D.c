/*******************************************************************
 *  2D HPC fit
 *  2023 AG HUth
 *  based on Uwe Kirst-> daytest, Oliver Radomski-> C++ implementation + 2D-histograms and Tobias Huth-> 2D fit
 *  FAU Erlangen-N�rnberg
 *******************************************************************/

#include <iostream>

#include "MultiLevel2D.h"

TMultiLevel2D::TMultiLevel2D(int starting_level, int level_increment, int number_of_levels)
{
  m_starting_level = starting_level;
  m_level_increment = level_increment;
  m_number_of_levels = number_of_levels;
}

void TMultiLevel2D::generate_matrix_2D_exp()
{

  int half_levels = int(m_number_of_levels / 2);
  double original_base = Daten.a_fit.i_null;   // baseline
  double original_amp = Daten.a_fit.i_channel; // amplitude

  // std::cout << "experimental: " << std::endl;
  // std::cout << "Daten.a_fit.i_null: " << Daten.a_fit.i_null << std::endl;
  // std::cout << "Daten.a_fit.i_channel: " << Daten.a_fit.i_channel << std::endl;
  // std::cout << "Daten.a_fit.sigma: " << Daten.a_fit.sigma << std::endl;
  double levels[half_levels];
  int filters_helper{int((half_levels - 1) / 2)};
  double range{0.8};

  // levels[0] = range;
  levels[int((half_levels - 1) / 2)] = 0;
  // levels[half_levels - 1] = -1 * range;

  for (int k = 0; k < filters_helper; k++)
  {
    levels[k] = range - k * (2.0 * range / (half_levels - 1));
    levels[half_levels - k - 1] = -1 * levels[k];
    // std::cout << "k: " << k << std::endl;

    // std::cout << levels[int(half_levels - 1 / 2) - k + 2] << std::endl;
  }

  // for (int j = 0; j < half_levels; j++)
  // {
  //   std::cout << "levels[i]: " << j << " " << levels[j] << std::endl;
  // }

  for (int i = 0; i < m_number_of_levels; i++) // m_number_of_levels levels
  {
    Daten.Dwell_2d_ptr = &Daten.Dwell_2d_MA[i];
    Daten.Dwell_2d_ptr->init((char *)"measured", Daten.Fitparameter.log_min_close,
                             Daten.Fitparameter.log_max_close,
                             Daten.Fitparameter.log_min_open,
                             Daten.Fitparameter.log_max_open,
                             Daten.Fitparameter.bins_per_log);

    if (i < m_number_of_levels / 2)
    {
      /***if (i == 0) {
        Daten.Dwell_2d_ptr = & Daten.Dwell_2d_MA[i];
        Daten.Dwell_2d_ptr -> init((char * )
          "measured", Daten.Fitparameter.log_min_close,
          Daten.Fitparameter.log_max_close,
          Daten.Fitparameter.log_min_open,
          Daten.Fitparameter.log_max_open,
          Daten.Fitparameter.bins_per_log);
      }***/
      Daten.a_fit.i_null = original_base + original_amp * levels[i];
      Daten.a_fit.i_channel = original_amp - original_amp * levels[i];

      // std::cout << "level_lower_rec: " << i << " : " << levels[i] << std::endl;
    }
    else
    {

      /***if (i == m_number_of_levels / 2) {
          Daten.Dwell_2d_ptr = & Daten.Dwell_2d_MA[i];
          Daten.Dwell_2d_ptr -> init((char * )
            "measured", Daten.Fitparameter.log_min_close,
            Daten.Fitparameter.log_max_close,
            Daten.Fitparameter.log_min_open,
            Daten.Fitparameter.log_max_open,
            Daten.Fitparameter.bins_per_log);

        }***/
      Daten.a_fit.i_null = original_base;
      Daten.a_fit.i_channel = original_amp - original_amp * levels[i - m_number_of_levels / 2];
      // std::cout << "level_upper_rec: " << levels[i - m_number_of_levels / 2] << std::endl;
    }
    // std::cout << "experimental: " << std::endl;
    Daten.hinkley(_HOHD);
    // std::cout<<"here2: "<<std::endl;
    Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);
    // std::cout<<"here3: "<<std::endl;
  }

  Daten.a_fit.i_null = original_base;   // return baseline to original value
  Daten.a_fit.i_channel = original_amp; // return amplitude to original value
}

double TMultiLevel2D::return_2d_matrix_LLH()
{

  int half_levels = int(m_number_of_levels / 2);
  double likelihood_2d{};
  double range{0.8};
  // std::cout<<"original: "<<original_sigma<<std::endl;

  double original_base = Daten.a_fit.i_null;
  double original_amp = Daten.a_fit.i_channel;

  // std::cout << "simulated: " << std::endl;
  // std::cout << "Daten.a_fit.i_null: " << Daten.a_fit.i_null << std::endl;
  // std::cout << "Daten.a_fit.i_channel: " << Daten.a_fit.i_channel << std::endl;
  // std::cout << "Daten.a_fit.sigma: " << Daten.a_fit.sigma << std::endl;
  double levels[half_levels];
  int filters_helper{int((half_levels - 1) / 2)};
  // std::cout<<"here4: "<<std::endl;

  levels[int((half_levels - 1) / 2)] = 0;
  // levels[half_levels - 1] = -1 * range;

  for (int k = 0; k < filters_helper; k++)
  {
    levels[k] = range - k * (2.0 * range / (half_levels - 1));
    levels[half_levels - k - 1] = -1 * levels[k];
    // std::cout << "k: " << k << std::endl;

    // std::cout << levels[int(half_levels - 1 / 2) - k + 2] << std::endl;
  }

  // for (int j = 0; j < half_levels; j++)
  // {
  //   std::cout << "levels[i]: " << j << " " << levels[j] << std::endl;
  // }

  for (int i = 0; i < m_number_of_levels; i++)
  { // m_number_of_levels levels

    Daten.Dwell_2d_ptr = &Daten.Dwell_2d_B;
    Daten.Dwell_2d_ptr->init((char *)"simulated", Daten.Fitparameter.log_min_close,
                             Daten.Fitparameter.log_max_close,
                             Daten.Fitparameter.log_min_open,
                             Daten.Fitparameter.log_max_open,
                             Daten.Fitparameter.bins_per_log);
    if (i < m_number_of_levels / 2)
    {
      /*** if (i == 0) {
         Daten.Dwell_2d_ptr = & Daten.Dwell_2d_B;
         Daten.Dwell_2d_ptr -> init((char * )
           "simulated", Daten.Fitparameter.log_min_close,
           Daten.Fitparameter.log_max_close,
           Daten.Fitparameter.log_min_open,
           Daten.Fitparameter.log_max_open,
           Daten.Fitparameter.bins_per_log);
       }***/

      Daten.a_fit.i_null = original_base + original_amp * levels[i];
      Daten.a_fit.i_channel = original_amp - original_amp * levels[i];

      // std::cout<<"level_lower: "<<levels[i]<<std::endl;
    }
    else
    {

      /***if (i == m_number_of_levels / 2) {
        Daten.Dwell_2d_ptr = & Daten.Dwell_2d_B;
        Daten.Dwell_2d_ptr -> init((char * )
          "simulated", Daten.Fitparameter.log_min_close,
          Daten.Fitparameter.log_max_close,
          Daten.Fitparameter.log_min_open,
          Daten.Fitparameter.log_max_open,
          Daten.Fitparameter.bins_per_log);

      }***/
      Daten.a_fit.i_null = original_base;
      Daten.a_fit.i_channel = original_amp - original_amp * levels[i - m_number_of_levels / 2];
      // std::cout<<"level_upper: "<<levels[i - m_number_of_levels/2]<<std::endl;
    }
    // std::cout << "simulated: " << std::endl;
    Daten.hinkley(_HOHD);
    Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);
    likelihood_2d += Daten.Dwell_2d_MA[i].lnlikelihood(Daten.Dwell_2d_B);
  }
  Daten.a_fit.i_null = original_base;   // return baseline to original value
  Daten.a_fit.i_channel = original_amp; // return amplitude to original value
  return likelihood_2d;
}
