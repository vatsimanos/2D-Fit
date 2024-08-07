#include "fft.h"
#include <algorithm>

namespace fft
{
    namespace helper
    {

        static double constexpr pi = 3.14159265358979323846;
        /*
        * ceiling(log2(x))
        * 
        * @return ceiling(log2(x))
        * 
        * @param x: log argument 
        */
        int upper_log2(const int x)
        {
            for (int i = 0; i < 30; ++i)
            {
                const int temp = 1 << i;
                if (temp >= x)
                {
                    return i;
                }
            }
            return 30;
        }

        /*
        *
        * generate the phase vector to assist fft
        * 
        * @return a vector of complex numbers that are uniformly located on the unity circle
        * 
        * @param len: is the number of points and length of the output vector. Must be 2^N
        */
        const std::vector<std::complex<double>> phase_vec(const int len)
        {
            std::vector<std::complex<double>> res(len);
            const double radius = 1;
            for (int i = 0; i < len; ++i)
            {
                const double phase = -2 * pi * i / len;
                res[i] = std::polar(radius, phase);
            }
            return res;
        }

        /*
        * butterfly forwarding the iteration of fft, Cooley Tukey algorithm
        * 
        * @param prev: fft in the previous iteration. 
        *               Is the time domain signal in the first iteration iff fft or frequency domain spectrum iff ifft
        * @param temp: fft in this iteration. 
        *               Is the frequency domina spectrum in the last iteration iff fft or time domain signal if ifft
        * @param phases: a vector of complex numbers that are uniformly located on the unity circle
        *  
        * @param turn: the iteration indicator, should be 0 if called outside
        * 
        * @param n_bits: log2(length_of_vectors), total number of iterations
        */
        void forward(std::vector<std::complex<double>> &prev, std::vector<std::complex<double>> &temp,
                     const std::vector<std::complex<double>> &phases, const int turn, const int n_bits)
        {
            if (turn == n_bits)
            {
                return;
            }

            const int group_size = 1 << (turn + 1); // size of butterfly group
            const int num_groups = prev.size() / group_size;
            const int phase_angular_freq = num_groups;
            for (int i_group = 0; i_group < num_groups; ++i_group)
            {
                const int base_index = i_group * group_size;
                // iterate through within the butterfly group
                for (int j = 0; j < group_size / 2; ++j)
                {
                    const int x0_index = base_index + j;
                    const int x1_index = base_index + group_size / 2 + j;
                    prev[x1_index] *= phases[j * phase_angular_freq];
                    temp[x0_index] = prev[x0_index] + prev[x1_index];
                    temp[x1_index] = prev[x0_index] - prev[x1_index];
                }
            }

            forward(temp, prev, phases, turn + 1, n_bits);
        }

        /*
        * A very fast O(N) bit reversal permutation algorithm. Published by John Wiley & Sons Ltd. 2002. 
        * https://www.researchgate.net/publication/227682626_A_new_superfast_bit_reversal_algorithm
        * 
        * @param vec: is the vector to permute. Must in the length of 2^N
        * @param n_bits: is N
        */
        void bit_reversal_permutation(std::vector<std::complex<double>> &vec, const int n_bits)
        {
            if (vec.size() <= 2)
            {
                return;
            }
            if (vec.size() == 4)
            {
                std::swap(vec[1], vec[3]);
                return;
            }
            std::vector<int> bit_rerversal(vec.size());
            // initialize the first 4 elements
            bit_rerversal[0] = 0;
            bit_rerversal[1] = 1 << (n_bits - 1);
            bit_rerversal[2] = 1 << (n_bits - 2);
            bit_rerversal[3] = bit_rerversal[1] + bit_rerversal[2];

            // theorm: for all nk = 2^k - 1 where k <= n_bits, we have
            //               bit_rerversal[nk] = bit_rerversal[n(k-1)] + 2^(n_bits-k)
            //               bit_rerversal[nk - i] ==  bit_rerversal[nk] - bit_rerversal[i]
            for (int k = 3; k <= n_bits; ++k)
            {
                const int nk = (1 << k) - 1;          // n_k = 2^k - 1
                const int n_km1 = (1 << (k - 1)) - 1; // n_(k-1) = 2^(k-1) - 1
                bit_rerversal[nk] = bit_rerversal[n_km1] + (1 << (n_bits - k));
                for (int i = 1; i <= n_km1; ++i)
                {
                    bit_rerversal[nk - i] = bit_rerversal[nk] - bit_rerversal[i];
                }
            }

            // permute the input vector according to bit_rerversal[]
            for (int i = 0; i < vec.size(); ++i)
            {
                if (bit_rerversal[i] > i)
                {
                    std::swap(vec[i], vec[bit_rerversal[i]]);
                }
            }
        }

    } // namespace helper

    std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &inputs)
    {
        if (inputs.empty())
        {
            return {};
        }
        const int n_bits = helper::upper_log2(inputs.size());
        const int len = 1 << n_bits;
        const std::vector<std::complex<double>> phases = helper::phase_vec(len);
        std::vector<std::complex<double>> prev(len);
        std::vector<std::complex<double>> temp(len);
        std::copy(inputs.begin(), inputs.end(), prev.begin());
        helper::bit_reversal_permutation(prev, n_bits);
        // butterfly forwarding from input to output, Cooley Tukey algorithm
        helper::forward(prev, temp, phases, 0, n_bits);
        return (n_bits % 2 == 1) ? temp : prev;
    }
    std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>> &inputs)
    {
        // How to calculate ifft by fft?
        // let x[n] denote time domain signal and X[k] denote frequency domain signal, then
        // x[n] = 1/N * sum(k=0..N-1) (X[k] * exp(1j*2*pi/N*k*n))

        // lets denote m = -k, then
        // x[n] = 1/N * sum(m=0..1-N) (X[m] * exp(-1j*2*pi/N*k*n)) == fft(X[m])
        
        // we know fft is circularly periodic, hence X[m] = X[-k] = X[N-k],
        // therefore we can flip the order of X[k] to get X[m]


        // flip the order of frequency spectrum
        std::vector<std::complex<double>> reverse_freq_spectrum(inputs);
        std::reverse(std::next(reverse_freq_spectrum.begin()), reverse_freq_spectrum.end());

        // normalization by multiplying 1/N to each element
        const double len = reverse_freq_spectrum.size();
        std::transform(reverse_freq_spectrum.begin(), reverse_freq_spectrum.end(), reverse_freq_spectrum.begin(),
                       [len](const std::complex<double> &num) { return num / len; });
        // fft
        return fft(reverse_freq_spectrum);
    }
    std::vector<int> round(const std::vector<std::complex<double>> &vec)
    {
        std::vector<int> res(vec.size());
        std::transform(vec.begin(), vec.end(), res.begin(), [](const std::complex<double> &num) -> int { return std::round(num.real()); });
        return res;
    }
    std::vector<double> real(const std::vector<std::complex<double>> &vec)
    {
        std::vector<double> res(vec.size());
        std::transform(vec.begin(), vec.end(), res.begin(), [](const std::complex<double> &num) -> double { return num.real(); });
        return res;
    }
} //namespace fft