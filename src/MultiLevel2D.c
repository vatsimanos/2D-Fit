/*******************************************************************
 *  2D HPC fit
 *  2023 AG HUth
 *  based on Uwe Kirst-> daytest, Oliver Radomski-> C++ implementation + 2D-histograms and Tobias Huth-> 2D fit
 *  FAU Erlangen-N�rnberg
 *******************************************************************/

#include <iostream>

#include "MultiLevel2D.h"


TMultiLevel2D::TMultiLevel2D(int starting_level, int level_increment, int number_of_levels) {
  m_starting_level = starting_level;
  m_level_increment = level_increment;
  m_number_of_levels = number_of_levels;
}

void TMultiLevel2D::generate_matrix_2D_exp() {
  int i;
  double original_base = Daten.a_fit.i_null;
  double original_amp = Daten.a_fit.i_channel;
  double levels[m_number_of_levels];

  for (int j = 0; j < int(m_number_of_levels / 2); j++) {
    levels[j] = j * (2.0 / m_number_of_levels) - 0.5;
  }

  for (i = 0; i < m_number_of_levels; i++) //m_number_of_levels levels
  {	
	        Daten.Dwell_2d_ptr = & Daten.Dwell_2d_MA[i];
        Daten.Dwell_2d_ptr -> init((char * )
          "measured", Daten.Fitparameter.log_min_close,
          Daten.Fitparameter.log_max_close,
          Daten.Fitparameter.log_min_open,
          Daten.Fitparameter.log_max_open,
          Daten.Fitparameter.bins_per_log);

    if (i < m_number_of_levels / 2) {
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
      //std::cout<<"level_lower_rec: "<<levels[i]<<std::endl;

    } else {

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
      //std::cout<<"level_upper_rec: "<<levels[i - m_number_of_levels/2]<<std::endl;
    }
	//std::cout<<"here1: "<<std::endl;
	Daten.hinkley(_HOHD);
	//std::cout<<"here2: "<<std::endl;
    Daten.Dwell_2d_ptr -> addjumps(Daten.Dwell_1d_ptr);
	//std::cout<<"here3: "<<std::endl;
  }
}

double TMultiLevel2D::return_2d_matrix_LLH() {
  int i;
  double likelihood_2d {};
  //std::cout<<"original: "<<original_sigma<<std::endl;

  double original_base = Daten.a_fit.i_null;
  double original_amp = Daten.a_fit.i_channel;
  double levels[m_number_of_levels];
	//std::cout<<"here4: "<<std::endl;
  for (int j = 0; j < int(m_number_of_levels / 2); j++) {
    levels[j] = j * (2.0 / m_number_of_levels) - 0.5;
  }

  for (i = 0; i < m_number_of_levels; i++) { //m_number_of_levels levels
  
	        Daten.Dwell_2d_ptr = & Daten.Dwell_2d_B;
        Daten.Dwell_2d_ptr -> init((char * )
          "simulated", Daten.Fitparameter.log_min_close,
          Daten.Fitparameter.log_max_close,
          Daten.Fitparameter.log_min_open,
          Daten.Fitparameter.log_max_open,
          Daten.Fitparameter.bins_per_log);
    if (i < m_number_of_levels / 2) {
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
      //std::cout<<"level_lower: "<<levels[i]<<std::endl;

    } else {
		
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
      //std::cout<<"level_upper: "<<levels[i - m_number_of_levels/2]<<std::endl;
    }
	Daten.hinkley(_HOHD);
	Daten.Dwell_2d_ptr -> addjumps(Daten.Dwell_1d_ptr);
	likelihood_2d += Daten.Dwell_2d_MA[i].lnlikelihood(Daten.Dwell_2d_B);
  }
  return likelihood_2d;
}