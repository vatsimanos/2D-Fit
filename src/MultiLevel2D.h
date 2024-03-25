/*******************************************************************
 *  2D HPC fit
 *  2023 AG HUth
 *  based on Uwe Kirst-> daytest, Oliver Radomski-> C++ implementation + 2D-histograms and Tobias Huth-> 2D fit
 *  FAU Erlangen-N�rnberg
 *******************************************************************/
#include "daten.h"


class TMultiLevel2D {

    private: int m_starting_level;
    int m_level_increment;
    int m_number_of_levels;

    public: TMultiLevel2D(int starting_level, int level_increment, int number_of_levels);

    void generate_matrix_2D_exp();

    double return_2d_matrix_LLH();

};