/*******************************************************************
 *  2D HPC fit
 *  2023 AG HUth
 *  based on Uwe Kirst-> daytest, Oliver Radomski-> C++ implementation + 2D-histograms and Tobias Huth-> 2D fit
 *  FAU Erlangen-Nürnberg
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
	int i;
		  
     	
  for (i=0;i < m_number_of_levels;i++)	
     {
     	Daten.Dwell_2d_ptr = &Daten.Dwell_2d_MA[i];
     	Daten.Dwell_2d_ptr->init((char *)"measured", Daten.Fitparameter.log_min_close,
        			Daten.Fitparameter.log_max_close, 
   	          Daten.Fitparameter.log_min_open, 
  				    Daten.Fitparameter.log_max_open, 
  				    Daten.Fitparameter.bins_per_log);    
  		Daten.a_fit.i_channel = m_starting_level + i*m_level_increment;		    
     	Daten.hinkley(_HOHD);		
     	Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);
     	//std::cout<<"level: "<<Daten.a_fit.i_channel<<" nr_of events: "<<Daten.Dwell_2d_ptr->events()<<std::endl;
     }		
}	
    
    
double TMultiLevel2D::return_2d_matrix_LLH()
{
	 int i;
	 double likelihood_2d{};

	 for (i=0;i < m_number_of_levels;i++)	
	 	{
	 	  Daten.Dwell_2d_ptr = &Daten.Dwell_2d_B;
  		Daten.Dwell_2d_ptr->init((char *)"simulated", Daten.Fitparameter.log_min_close,
  				        Daten.Fitparameter.log_max_close, 
  				        Daten.Fitparameter.log_min_open, 
  				        Daten.Fitparameter.log_max_open, 
  				        Daten.Fitparameter.bins_per_log);
  	  Daten.a_fit.i_channel = m_starting_level + i*m_level_increment;		    
     	Daten.hinkley(_HOHD);		
     	Daten.Dwell_2d_ptr->addjumps(Daten.Dwell_1d_ptr);
     	likelihood_2d += Daten.Dwell_2d_MA[i].lnlikelihood(Daten.Dwell_2d_B);	   
     	//std::cout<<"level: "<< m_starting_level + i*m_level_increment<<" likelihood_2d: "<<likelihood_2d<<std::endl;   
	  }	
 return likelihood_2d;	  
}
