/*

*/

#ifndef SSPE_H_INCLUDED
#define SSPE_H_INCLUDED

#define _USE_MATH_DEFINES

#include <BasicJuceHeader.h>
#include <Eigen\Dense>
#include <math.h>
#include <random>
#include <complex>
#include <chrono>

#include <iostream>
#include <fstream>

#include <limits>
#include <iomanip>
#include <sstream>
typedef std::numeric_limits< double > dbl;

//using namespace StatePhaseEst;

namespace StatePhaseEst
{ 
    using vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;

    enum SSPE_PARAM
    {
        DATA_FS,
        STRIDE,
        OBS_ERR_EST_SSPE,
        Q_EST_ONE_SSPE,
        Q_EST_TWO_SSPE,
        Q_EST_THREE_SSPE
    };

    // generate variables to guess on quality of phase esimtation
    struct normal_random_variable
    {
        
        normal_random_variable(Eigen::MatrixXd const& covar)
            : normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
        {}

        normal_random_variable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
            : mean(mean)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
            transform = eigenSolver.eigenvectors() * (eigenSolver.eigenvalues().cwiseSqrt().asDiagonal());
        }

        Eigen::VectorXd mean;
        Eigen::MatrixXd transform;

        Eigen::VectorXd operator()() const
        {
            static std::mt19937 gen{ std::random_device{}() };
            static std::normal_distribution<> dist;

            return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
            //return Eigen::VectorXd();
        }
        
    };

    class SSPE {
    public:
        SSPE() 
        :   fs          (1000)
        ,   foiIndex    (0)
        ,   sigmaObs    (0)
        ,   stride      (30000/fs)
        ,   histSize(0)
        ,   dataFs(30000)
        ,   ampEst(vector(1))
        ,   qEst(vector(1))
        ,   obsErrorEst(1)
       // ,   myfile(std::ofstream("D:\\TNEL\\oep-installation\\state-phase-est\\StatePhaseEst\\Source\\bufdat.txt", std::ios::app))

        {
            ampEst(0) = 0.999995434106301;
            qEst(0) = 50;

            allX = Array<Matrix>();
            allP = Array<Matrix>();
            dsBuf = Array<float>();
            prevBuf = Array<float>();
            complexArray = Array<std::complex<double>>();
        }
        
        ~SSPE() { }

        void reset()
        {
            hasBeenUsed = false;
            //freqsSet = false;
        }

        bool hasBeenFit()
        {
            return hasBeenUsed;
        }

        // returns true if successful.
        bool setParams(int paramIndex, double value)
        {
            switch (paramIndex)
            {
            case DATA_FS:
                dataFs = value;
                //stride = dataFs / fs;
                break;
            case STRIDE:
                stride = value;
                break;
            case OBS_ERR_EST_SSPE:
                obsErrorEst = value;
                break;
            case Q_EST_ONE_SSPE:
                qEst(0) = value;
                break;
            case Q_EST_TWO_SSPE:
                qEst(1) = value;
                break;
            case Q_EST_THREE_SSPE:
                qEst(2) = value;
                break;
            }
           
            return true;
        }

        bool setHistSize(int histSize)
        {
            histSize = histSize;

            allPEB = Array<Matrix>();
            allXEB = Array<Matrix>();

            for (int i = 0; i <= histSize; i+=stride)
            {
                allXEB.add(Matrix());
                allPEB.add(Matrix());
            }

            return true;
        }

        bool setFreqs(Array<float> newFreqs)
        {
            freqs = vector(newFreqs.size());
            for (int i = 0; i < newFreqs.size(); i++)
            {
                freqs(i) = newFreqs[i];
            }
            return true;
        }

        const int getFs()
        {
            return fs;
        }

        bool setDesFreqIndex(int foi)
        {
            int nFreqs = freqs.size();
         
            if (foi <= nFreqs)
            {
                foiIndex = foi;
                freqsSet = true;
                return freqsSet;
            }

            freqsSet = false;
            return freqsSet;
        }

        Array<std::complex<double>> evalBuffer(const float* rpBuf, int nSamples)
        { 

            int arraySize = ((nSamples-1) / stride) + 1; // nSamples -1 in case perfectly divisible by stride so we round down. +2 to account for overflow value

            // clear allP
            allX.clearQuick();
            if (allX.size() != arraySize+1) // Plus one because it also needs to hold prev state
            {
                allX.resize(arraySize+1);
            }
            allX.set(0,prevState); // Append prevState
            
            // clear allP
            allP.clearQuick(); 
            if (allP.size() != arraySize+1) // Plus one because it also needs to hold prev cov
            {
                allP.resize(arraySize+1);
            }
            allP.set(0,prevCov); // Append prevCov

           // clear dsBuf
           dsBuf.clearQuick();
           if (dsBuf.size() != arraySize)
           {
               dsBuf.resize(arraySize);
           }

           // Add samples to dsBuf from rpBuf and prevBuf (if needed for overflow)
            for (int i = nSamples, n = arraySize-1; n>=0; i-=stride, n--)
            {
                if (i < 0)
                {
                    dsBuf.set(n, prevBuf[nSamples + i - 1]); // go back i samples from the end (i will be negative because we passed the 0th sample of the current buffer
                    // and are going to the prev buffer.
                    std::cout << "should never get here" << std::endl;
                    std::cin.get();
                }
                else
                {
                    dsBuf.set(n, rpBuf[i - 1]);
                }
            }
                        
            // Start setting allX at index 1 because 0th carried over from previous buffer
            for (int i = 1; i <= arraySize; i++)
            {
                Matrix newState;
                Matrix newStateCov;
                oneStepKalmanEval(newState, newStateCov, allX[i - 1], dsBuf[(i - 1)], allP[i - 1], phi, Q, sigmaObs, M);
                allX.set(i,newState);
                allP.set(i,newStateCov);
            }
            
            complexArray.clearQuick();
            if (complexArray.size() != arraySize)
            {
                complexArray.resize(arraySize);
            }
            prevState = allX.getLast();
            prevCov = allP.getLast();

            for (int i = 1; i <= arraySize; i++)
            {
                 complexArray.set(i-1,std::complex<double>(allX[i](foiIndex, 0), allX[i](foiIndex + 1, 0)));
            }

            // set prevBuf (maybe change to only save latest)
            prevBuf.clearQuick();
            if (prevBuf.size() != nSamples)
            {
                prevBuf.resize(nSamples);
            }
            for (int i = 0; i < nSamples; i++)
            {
                prevBuf.set(i, rpBuf[i]);
            }

            return complexArray;      
        }

        void oneStepKalmanEval(Matrix& newState, Matrix& newStateCov, Matrix state, float measuredVal, Matrix stateCov, Matrix phi, Matrix Q, float sigmaObs, Matrix M)
        {
            Matrix stateEst = phi * state;
            Matrix pEst = phi * stateCov * phi.transpose() + Q;

            Matrix kalmanGain = (pEst * M) / ((M.transpose() * pEst * M)(0) + sigmaObs);

            Matrix temp = M.transpose() * stateEst;
            newState = stateEst + kalmanGain * (measuredVal - temp(0));
            newStateCov = pEst - kalmanGain * M.transpose() * pEst;
        }


        float prctile(float percent, Array<float> arr)
        {
            
            int nArr = arr.size();
            int pctIndex = std::ceil(nArr * percent / 100);
            return arr[pctIndex];
            
            return 0;
        }

        void genKalmanParams(Matrix& phi, Matrix& Q, Matrix& M, vector freqs, vector ampVec, vector sigmaFreqs)
        {
            int nFreqs = freqs.size();
            Array<Matrix> rotMat = Array<Matrix>(); // R(w)
            Array<Matrix> varMat = Array<Matrix>(); // Q

            for (int i = 0; i < nFreqs; i++)
            {
                // Add rotation matrix
                double freq = freqs(i);
                Matrix mat = Matrix(2, 2);
                double pi = 3.141592653589793;
                mat(0, 0) = cos(2 * pi * freq / fs);
                mat(0, 1) = -1 * sin(2 * pi * freq / fs);
                mat(1, 0) = sin(2 * pi * freq / fs);
                mat(1, 1) = cos(2 * pi * freq / fs);
                Matrix ampMat = mat * ampVec(i);
                rotMat.add(ampMat);

                // Add variance matrix
                mat = Matrix(2, 2);
                mat(0, 0) = sigmaFreqs(i);
                mat(0, 1) = 0;
                mat(1, 0) = 0;
                mat(1, 1) = sigmaFreqs(i);
                varMat.add(mat);
            }

            phi = Matrix::Zero(nFreqs * 2, nFreqs * 2);
            Q = Matrix::Zero(nFreqs * 2, nFreqs * 2);

            // Add rotMat and varMat to phi and Q along the diagonal to improve speed of processing
            for (int i = 0; i < nFreqs; i++)
            {
                // Zeroes before mat
                for (int j = 0; j < i * 2; j++)
                {
                    //phi(i, j) = 0;
                    //phi(i + 1, j) = 0;
                    //Q(i, j) = 0;
                    //Q(i + 1, j) = 0;
                }

                // Put mat along diagonal
                phi(i * 2, i * 2) = rotMat[i](0, 0);
                phi(i * 2, i * 2 + 1) = rotMat[i](0, 1);
                phi(i * 2 + 1, i * 2) = rotMat[i](1, 0);
                phi(i * 2 + 1, i * 2 + 1) = rotMat[i](1, 1);

                Q(i * 2, i * 2) = varMat[i](0, 0);
                Q(i * 2, i * 2 + 1) = varMat[i](0, 1);
                Q(i * 2 + 1, i * 2) = varMat[i](1, 0);
                Q(i * 2 + 1, i * 2 + 1) = varMat[i](1, 1);

                // Pad zeroes after mat
                for (int j = i * 2 + 2; j < nFreqs * 2; j++)
                {
                    //phi(i * 2, j) = 0;
                    //phi(i * 2 + 1, j) = 0;
                    //Q(i * 2, j) = 0;
                    //Q(i * 2 + 1, j) = 0;
                }
            }
            
            //M = reshape([ones(length(freqs), 1), zeros(length(freqs), 1)]', length(freqs)*2,1); % for extracting real part of signal from analytic 'x'
            M = Matrix::Zero(nFreqs * 2,1);
            for (int i = 0; i < nFreqs * 2; i++)
            {
                if (i % 2 == 0)
                {
                    M(i,0) = 1;
                }
                //else
                //{
                    //M(i,0) = 0.0;
                //}
            }

        }

        void fitModel(Array<double> data_array)
        {
            std::cout.precision(dbl::max_digits10);
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            //int stride = 40000 / fs;
            int stride = dataFs / fs;
            std::ofstream myfile("E:\\TNEL\\oep-installation\\state-phase-est\\StatePhaseEst\\Source\\bufdat_justsin_3101024.csv", std::ios::app);
            //std::ofstream myfile("C:\\Users\\Ephys\\Documents\\Github\\OE\\oecmake\\phase-calculator\\StatePhaseEst\\Source\\bufdat_justsin_3101024.csv", std::ios::app);
           // std::cout << "Stride: " << stride << std::endl;
            //std::cout << "data Size: " << data_array.size() << std::endl;
            const double* inputSeries = data_array.begin();
            vector data = vector(data_array.size() / stride);
            for (int i = 0; i < data_array.size()/stride; i++)
            {
                //std::stringstream ss;
                //ss << std::setprecision( std::numeric_limits<double>::digits10+2);
                //ss << std::setprecision(std::numeric_limits<double>::max());
                //ss << inputSeries[i * stride];
                myfile << std::fixed << std::setprecision(dbl::max_digits10 + 2) << inputSeries[i * stride] << ",";
                data(i) =  inputSeries[i*stride];
                if (i < 3)
                {
                    std::cout << i << " : " << data(i) << std::endl;
                }      
            }
            
            myfile << "\n" << std::endl;
            //myfile.close();
            //std::cin.ignore();
            if (freqsSet == false)
            {
                // Freqs have not been set yet. need to set before moving on.
                jassertfalse;
                return;
            }
            
            int nSamples = data.size();
            int nFreqs = freqs.size();
            Matrix xstart = Matrix::Zero(2 * nFreqs, 1);
            Matrix Pi;
            Pi.setIdentity(2 * nFreqs, 2 * nFreqs);
            Matrix P = 0.001 * Pi;

            Matrix tempPhi;
            Matrix tempQ;
            Matrix tempM;
            double tempSigmaObs;
            vector ampEst;
            vector allQ;
            vector freqEst;
            vector omega;
            
            // Create initial phi, Q and M based on params used last ( or estimates if first run)
            if (hasBeenFit())
            {
                tempPhi = phi;
                tempQ = Q;
                tempM = M;
                tempSigmaObs = sigmaObs;
                ampEst = ampVec;
                allQ = sigmaFreqs;
                freqEst = freqs / fs;
            }
            else
            {
                // Run AR model

                // Just hardcode for now?
                //initParams(init_f, init_a, init_sigma, init_sigma_obs);
                /*
                ampEst = vector(1); // = // FROM MATLAB (0.9991)
                ampEst(0) = 0.999995434106301;
                allQ = vector(1); // = // from matlab (30.8456)
                allQ(0) = 50;
                tempSigmaObs = obsErrorEst; // = // from matlab
                freqEst = vector(1);
                freqEst(0) = 5 / fs;
                freqs(0) = 5;
                
                
                ampEst = vector(3); // = // FROM MATLAB (0.9991)
                ampEst(0) = 0.99;
                ampEst(1) = 0.99;
                ampEst(2) = 0.99;
                allQ = vector(3); // = // from matlab (30.8456)
                allQ(0) = 50;
                allQ(1) = 50;
                allQ(2) = 50;
                tempSigmaObs = 1; // = // from matlab
                
                freqs = vector(3);
                freqs(0) = .3;
                freqs(1) = 10;
                freqs(2) = 60;
                
                freqEst = vector(3);
                freqEst = freqs / fs;
                */

                tempSigmaObs = obsErrorEst;

                nFreqs = freqs.size();
                ampEst.resize(nFreqs);
                allQ.resize(nFreqs);
                freqEst.resize(nFreqs);


                for (int i = 0; i < nFreqs; i++)
                {
                    ampEst(i) = 0.99;
                    allQ(i) = 50;
                    freqEst(i) = freqs[i] / fs;
                }
                Matrix Pi;
                Pi.setIdentity(2 * nFreqs, 2 * nFreqs);
                prevCov = 0.001 * Pi;
                P = 0.001 * Pi;
                //prevCov = Matrix::Zero(2*freqs.size(), 2*freqs.size()); // init prevCov and prevState to 0s
                prevState = Matrix::Zero(2 * nFreqs, 1);
                //tempQ = P;
                xstart = Matrix::Zero(2 * nFreqs, 1);
                
            }

            //vector freqEst = freqs / fs;
            
            genKalmanParams(tempPhi, tempQ, tempM, freqs, ampEst, allQ); // gen initial kalman params
            

            int iter = 1;
            float errorVal = FLT_MAX;
            Matrix prevStateModel = prevState;
            Matrix prevCovModel = prevCov;

 

            Array<Matrix> allP_N_N1 = Array<Matrix>();

            Array<Matrix> allXSmooth = Array<Matrix>();
            Array<Matrix> allPSmooth = Array<Matrix>();
            Array<Matrix> allJSmooth = Array<Matrix>();

            while (iter < 400 && errorVal > 1e-3)
            {
                for (int i = 0; i < nSamples; i++)
                {
                    allP_N_N1.add(Matrix::Zero(2 * nFreqs, 2 * nFreqs));
                }

                allXSmooth = Array<Matrix>();
                allPSmooth = Array<Matrix>();
                allJSmooth = Array<Matrix>();
                for (int i = 0; i < nSamples; i++)
                {
                    allXSmooth.add(Matrix::Zero(2 * nFreqs, 1));
                    allPSmooth.add(Matrix::Zero(2 * nFreqs, 2 * nFreqs));
                    if (i != nSamples - 1)
                    {
                        allJSmooth.add(Matrix::Zero(2 * nFreqs, 2 * nFreqs));
                    }
                }
                // run forward kalman filter
                Array<Matrix> allX = Array<Matrix>();
                allX.add(prevStateModel); // Append prevState

                Array<Matrix> allP = Array<Matrix>();
                allP.add(prevCovModel); // Append prevCov

                //std::cout << "tmpsgmaobs: " << tempSigmaObs << std::endl;
                //std::cout << "tempPhi b4:\n " << tempPhi << std::endl;
                //std::cout << "tempQ b4:\n " << tempQ << std::endl;
                //std::cout << " tempM b4:\n " << tempM << std::endl;

                //std::cout << "init x: " << allX[0] << std::endl;
                //std::cout << "init P: " << allP[0] << std::endl;
                //std::cout << "nsamp: " << nSamples << std::endl;
                //std::cin.ignore();
                for (int i = 0; i < nSamples; i++)
                {
                    Matrix newState;
                    Matrix newStateCov;
                    oneStepKalman(newState, newStateCov, allX[i], data(i), allP[i], tempPhi, tempQ, tempSigmaObs, tempM);
                    allX.add(newState);
                    allP.add(newStateCov);
                    double real = allX[i](0);
                    double img = allX[i](1);
                    //std::cout << i << " : " << newStateCov << std::endl;

                }
                //std::cin.ignore();
                allX.remove(0); // remove prevState
                allP.remove(0); // remove prevCov

                // Run fixed interval smoother - Shumway, Stoffer, 82
                allXSmooth.set(nSamples-1, allX[nSamples-1]);
                allPSmooth.set(nSamples-1, allP[nSamples-1]);
               // allJSmooth = Array<Matrix>();
                //allXSmooth= allX;
                //allPSmooth = allP;
               // allJSmooth = allP;
                //allJSmooth.remove(nSamples - 1);

                for (int i = nSamples - 1; i >= 0; i--)
                {
                    Matrix smoothState;
                    Matrix smoothCov;
                    Matrix smoothJ;
                    int nxs = allX.size();
                    int nps = allP.size();
                    if (i == nSamples - 1)
                    {
                        smoother(smoothState, smoothCov, smoothJ, allXSmooth[i], allX[i], allPSmooth[i], allP[i], tempPhi, tempQ);
                        //std::cout << i << " :allX[i]: " << allX[i] << std::endl;
                    }
                    else
                    {
                        smoother(smoothState, smoothCov, smoothJ, allXSmooth[i + 1], allX[i], allPSmooth[i + 1], allP[i], tempPhi, tempQ);
                        //std::cout << i << " :allX[i+`1]: " << allX[i+1] << std::endl;
                    }
                    //std::cout << i << ", smoothState: " << smoothCov << std::endl;
                    
                    //std::cout << "allXSmooth[i + 1]: " << allXSmooth[i + 1] << std::endl;
                   
                    //allXSmooth[i] = smoothState;
                    //allPSmooth[i] = smoothCov;
                    //allJSmooth[i] = smoothJ;
                    allXSmooth.set(i, smoothState);
                    allPSmooth.set(i, smoothCov);
                    allJSmooth.set(i, smoothJ);

                }
                //std::cin.ignore();
                /*
                std::cout << "allXSmooth: \n" << allXSmooth[0] << std::endl;
                std::cout << "allXSmooth end: \n" << allXSmooth[nSamples - 1] << std::endl;
                std::cout << "allPSmooth: \n" << allPSmooth[0] << std::endl;
                std::cout << "allPSmooth end: \n" << allPSmooth[nSamples-1] << std::endl;
                std::cout << "allJSmooth: \n" << allJSmooth[0] << std::endl;
                std::cout << "allJSmooth end: \n" << allJSmooth[nSamples - 2] << std::endl;
                std::cin.get();
                */
                // need to estimate P_t_(t-1) for t = 1:N
                Matrix PTemp = tempPhi * allP[nSamples-1] * tempPhi.transpose() + tempQ;
                double PArrayForm = (tempM.transpose() * PTemp * tempM)(0);
                double PTransformedMatrix = 1.0 / (PArrayForm + tempSigmaObs);
                Matrix K = PTemp * tempM * PTransformedMatrix;
                Matrix Identity = Matrix().Identity(PTemp.rows(), PTemp.cols());
                Matrix P_N_N1 = (Identity - K * tempM.transpose()) * tempPhi * allP[nSamples - 1];
                
                allP_N_N1.set(nSamples-1, P_N_N1);
                //std::cout << nSamples-1 << ", all_PNN1: " << P_N_N1 << std::endl;

                for (int i = nSamples - 2; i > 0; i--)
                {
                    if (i == nSamples - 1)
                    {
                        P_N_N1 = allP[i] * allJSmooth[i - 1].transpose() + allJSmooth[i] * (allP_N_N1[i] - tempPhi * allP[i]) * allJSmooth[i - 1].transpose();
                    }
                    else
                    {
                        P_N_N1 = allP[i] * allJSmooth[i - 1].transpose() + allJSmooth[i] * (allP_N_N1[i + 1] - tempPhi * allP[i]) * allJSmooth[i - 1].transpose();
                    }
                    
                    allP_N_N1.set(i, P_N_N1);
                    //std::cout << i << ", all_PNN1: " << P_N_N1 << std::endl;
                }
                Matrix eyee;
                allP_N_N1.set(0, eyee.setIdentity(2 * nFreqs, 2 * nFreqs));
                //std::cout << "0" << ", all_PNN1: " << allP_N_N1[0] << std::endl;
                //std::cout << "allpnn1: " << allP_N_N1[0] << std::endl;
                //std::cout << "allpnn1 end: " << allP_N_N1[nSamples-1] << std::endl;
                //std::cin.ignore();

                // EM Alg from Soulat, 2019
                Matrix A = P + xstart * xstart.transpose();
                //std::cout << "pre"<< ": " << A << std::endl;
                Matrix B = Matrix::Zero(allPSmooth[0].rows(), allPSmooth[0].cols());
                Matrix C = Matrix::Zero(allPSmooth[0].rows(), allPSmooth[0].cols());
                tempSigmaObs = 0;

                for (int i = 0; i < nSamples; i++)
                {
                    if (i > 0)
                    {
                        A = A + allPSmooth[i - 1] + allXSmooth[i - 1] * allXSmooth[i - 1].transpose();
                        //std::cout << i << ": " << A << std::endl;
                        //std::cout << i << "allPSmooth[i - 1]: " << allPSmooth[i - 1] << std::endl;
                        //std::cout << i << "allXSmooth[i - 1]: " << allXSmooth[i - 1] << std::endl;
                        //std::cout << i << "allXSmooth[i - 1]  * allXSmooth[i - 1].transpose: " << allXSmooth[i - 1] * allXSmooth[i - 1].transpose() << std::endl;
                        //std::cout << i << "allPSmooth[i - 1] + allXSmooth[i - 1] * allXSmooth[i - 1].transpose(): " << allPSmooth[i - 1] + allXSmooth[i - 1] * allXSmooth[i - 1].transpose() << std::endl;
                        //std::cout << i << "allXSmooth[i - 1].transpose: " << allXSmooth[i - 1].transpose() << std::endl;
                        //std::cout << i << "allXSmooth[i - 1].transpose: " << allXSmooth[i - 1].transpose() << std::endl;
                        //std::cin.ignore();
                        //std::cout << "all_PNNN1 cols: " << allP_N_N1[i].cols() << "\nall_PNN1 rows: " << allP_N_N1[i].rows() << std::endl;
                        //std::cout << "allx cols: " << allXSmooth[i].cols() << "\nallx rows: " << allXSmooth[i].rows() << std::endl;
                        B = B + allP_N_N1[i] + allXSmooth[i] * allXSmooth[i - 1].transpose();
                    }
                    C = C + allPSmooth[i] + allXSmooth[i] * allXSmooth[i].transpose();

                    // MS - Yikes... So you can only add scalar values to Arrays for some reason.
                    // MS - Also it doesn't let you treat the product of two Matrixes as a double, even if it will be. so temp is a 1x1 matrix and we just pull out the single number
                    //R = R+ M'* squeeze(newAllP(:,:,i)) * M + (y(i) - M'*newAllX(:,i)) *(y(i) - M'*newAllX(:,i))' ;
                    double temp = (tempM.transpose() * allPSmooth[i] * tempM)(0);
                    double temp1 = (tempM.transpose() * allXSmooth[i])(0);
                    //double temp2 = (tempM.transpose() * allXSmooth[i])
                    double temp2 = (data(i) - temp1);
                    double temp3 = temp2 * temp2;
                    //Eigen::ArrayXXf temp2 = (tempM.transpose() * allXSmooth[i]).transpose();
                    //std::cout << "i: " << i << std::endl;
                    //std::cout << i << std::endl;
                    //std::cout << "tempSigmaObsb4: " << tempSigmaObs << std::endl;
                    if (i == 0)
                    {
                        //std::cout << "hello" << std::endl;
                        tempSigmaObs = temp + temp3;
                    }
                    else
                    {
                        tempSigmaObs = tempSigmaObs + temp + temp3;
                    }
                    
                    //std::cout << "tempM: " << tempM << std::endl;
                    //std::cout << "appPsmooth: " << allPSmooth[i] << std::endl;
                    //std::cout << "allXSmooth: " << allXSmooth[i] << std::endl;
                    //std::cout << "data(i): " << data(i) << std::endl;
                    //std::cout << "tempSigmaObs2: " << tempSigmaObs << std::endl;

                    //std::cout << "temp: " << temp << std::endl;
                    //std::cout << "temp1: " << temp1 << std::endl;
                    //std::cout << "temp2: " << temp2 << std::endl;
                    //std::cout << "temp3 : " << temp3 << std::endl;
                    //std::cout << "sampR: " << temp + temp3 << std::endl;
                    //std::cout << "tempSigmaObsafter: " << tempSigmaObs << std::endl;
                    //std::cin.ignore();
                }

                tempSigmaObs = (1.0 / double(nSamples)) * tempSigmaObs;
                //std::cout << "a: " << A << std::endl;
                //std::cout << "b: " << B << std::endl;
                //std::cout << "c: " << C << std::endl;
                //std::cout << "tempsobs : " << tempSigmaObs << std::endl;
                //std::cout << "nsam: " << double(nSamples) << std::endl;
               // std::cin.ignore();

                

                
                
                vector oldFreq = freqEst * 1000 / (2 * M_PI);

                freqEst = vector::Zero(nFreqs);
                ampEst = vector::Zero(nFreqs);
                allQ = vector::Zero(nFreqs);

                for (int i = 0; i < nFreqs; i++)
                {

                    Matrix ATmp = extract2x2(A, i);
                    Matrix BTmp = extract2x2(B, i);
                    Matrix CTmp = extract2x2(C, i);
                    
                    //std::cout << "after genK ATmp at end:\n " << ATmp << std::endl;
                    //std::cout << "after genK BTmp at end:\n " << BTmp << std::endl;
                    //std::cout << "(BTmp(1, 0) - BTmp(0, 1)) / BTmp.trace()):\n " << (BTmp(1, 0) - BTmp(0, 1)) / BTmp.trace() << std::endl;
                    //std::cout << "after genK CTmp at end:\n " << CTmp << std::endl;
                    //std::cin.ignore();

                    freqEst(i) = (atan((BTmp(1, 0) - BTmp(0, 1)) / BTmp.trace()));
                    ampEst(i) = sqrt(pow((BTmp(1, 0) - BTmp(0, 1)), 2) + pow(BTmp.trace(), 2)) / ATmp.trace();
                    //std::cout << "CTmp.trace(): " << CTmp.trace() << std::endl;
                    //std::cout << "pow(ampEst(i), 2): " << pow(ampEst(i), 2) << std::endl;
                    //std::cout << " ATmp.trace(): " << ATmp.trace() << std::endl;
                    //std::cout << " (float(2.0 * nSamples)): " << (float(2.0 * nSamples)) << std::endl;
                    //std::cout << "CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace(): " << CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace() << std::endl;
                    //std::cout << "(float(2.0 * nSamples)) *CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace(): " << (float(2.0 * nSamples)) * (CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace()) << std::endl;
                    allQ(i) = 1.0 / (double(2.0 * nSamples)) * (CTmp.trace() - pow(ampEst(i), 2) * ATmp.trace());

                    //// ALERT CHECK AMP EST POW FUNCTION. WE NEED TO DO .^2 vs other ones are just ^2
                    // Matlab functions below
                    //ampEst(numFreqs) = sqrt((B_tmp(2,1) - B_tmp(1,2))^2 + trace(B_tmp)^2)/trace(A_tmp);
                    //allQ(numFreqs) = 1/(2*length(y)) * (trace(C_tmp) - ampEst(numFreqs).^2 * trace(A_tmp));

                }
                /*
                std::cout << "oldfreq: " << oldFreq << std::endl;
                std::cout << "tmpsgmaobs: " << tempSigmaObs << std::endl;
                std::cout << "freqest: " << freqEst << std::endl;
                std::cout << "ampEst: " << ampEst << std::endl;
                std::cout << "allQ: " << allQ << std::endl;
                */
               // std::cin.ignore();
                

                omega = freqEst * 1000.0 / (2.0 * M_PI);

                genKalmanParams(tempPhi, tempQ, tempM, omega, ampEst, allQ);
                //std::cout << "tempQ:\n " << tempQ << std::endl;
                //std::cout << "after genK allQ at end:\n " << allQ << std::endl;
                //std::cin.get();


                /*
                std::cout << "tmpsgmaobs: " << tempSigmaObs << std::endl;
                std::cout << "after genK tempPhi at end:\n " << tempPhi << std::endl;
                std::cout << "after genK ampEst at end:\n " << ampEst << std::endl;
                std::cout << "after genK freqEst at end:\n " << freqEst << std::endl;
                std::cout << "after genK tempQ at end:\n " << tempQ << std::endl;
                std::cout << "after genK omega at end:\n " << omega << std::endl;
                std::cout << "----------------------------" << std::endl;
                */
                ///std::cin.ignore();
                
                // Reset state vector
                //allX.clearQuick();
                //allX.add(allXSmooth.getLast()); // Append prevState

                //allP.clearQuick();
                //allP.add(allPSmooth.getLast()); // Append prevCov

                iter += 1;
                //  MAKE SURE cwiseAbs does what we think it does
                //std::cout << "omega: "  << omega << std::endl;
                //std::cout << "oldFreq: "  << oldFreq << std::endl;
                errorVal = ((omega - oldFreq).cwiseAbs()).sum();
                //std::cout << "errorval: " << errorVal << std::endl;

                //std::cin.get();
            }

            // Got good values send them to main thread. Maybe do this in run function

            // getLock()
            // {
            phi = tempPhi;
            Q = tempQ;
            M = tempM;
            sigmaObs = tempSigmaObs;
            ampVec = ampEst;
            sigmaFreqs = allQ;
            freqs = omega;
            hasBeenUsed = true;
            //std::cout << "freqs:\n" << freqs << std::endl;

            //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            
            std::cout << "phi:\n" << phi << std::endl;
            std::cout << "Q:\n" << Q << std::endl;
            std::cout << "M:\n" << M << std::endl;
            std::cout << "sigmaObs:\n" << sigmaObs << std::endl;
            std::cout << "ampVev:\n" << ampVec << std::endl;
            std::cout << "sigmaFreqs:\n" << sigmaFreqs << std::endl;
            std::cout << "freqs:\n" << freqs << std::endl;
            std::cout << "iter: " << iter << std::endl;

            std::cout << "---------------------------" << std::endl;
            
            std::cout << "fit!!" << std::endl;
            //std::cin.get();
            
            // }

            
           // std::cout << "Time difference fit model = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

        }

        //Eigen::MatrixXf extract2x2(Matrix inputMatrix, int i)
        Matrix extract2x2(Matrix inputMatrix, int i)
        {
            
            int rowcol = i * 2;
            Matrix x = Matrix::Zero(2, 2);
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    x(j, k) = inputMatrix(rowcol + j, rowcol + k);
                }
            }

            return x;
            
            //return Matrix::Zero(2, 2);
        }


        void smoother(Matrix& smoothX, Matrix& smoothP, Matrix& smoothJ, Matrix x, Matrix x_prev, Matrix P, Matrix P_prev, Matrix phi, Matrix Q)
        {
            // Matlab code uses x and x_n, swapped to use x=x_prev and x_n=x. Ditto for P
            Matrix P_one = phi * P_prev * phi.transpose() + Q;

            smoothJ = P_prev * phi.transpose() * P_one.inverse();
            smoothX = x_prev + smoothJ * (x - phi * x_prev);
            smoothP = P_prev + smoothJ * (P - P_one) * smoothJ.transpose();

            
            //std::cout << "smoothX:\n" << smoothX << std::endl;
            /*
            std::cout << "smoothP:\n" << smoothP << std::endl;
            std::cout << "smoothJ:\n" << smoothJ << std::endl;
            std::cout << "x:\n" << x << std::endl;
            std::cout << "x_prev:\n" << x_prev << std::endl;
            std::cout << "P:\n" << P << std::endl;
            std::cout << "P_prev:\n" << P_prev << std::endl;
            std::cout << "phi:\n" << phi << std::endl;
            std::cout << "Q:\n" << Q << std::endl;
            std::cout << "P_one.inverse:\n" << P_one.inverse() << std::endl;
            std::cout << "P_one:\n" << P_one << std::endl;
            std::cin.get(); 
            */
            
        }

        void oneStepKalman(Matrix& newState, Matrix& newStateCov, Matrix state, float measuredVal, Matrix stateCov, Matrix phi, Matrix Q, double sigmaObs, Matrix M)
        {
            //float breakme = phi(13, 5364);
            Matrix stateEst = phi * state;
            Matrix pEst = phi * stateCov * phi.transpose() + Q;

            //int ncol = M.size();
            //int nlegn = phiEst.size();

            //float temp = (M.transpose() * pEst * M)(0) + sigmaObs;
            Matrix kalmanGain = (pEst * M) / ((M.transpose() * pEst * M)(0) + sigmaObs);

            Matrix temp = M.transpose() * stateEst;
            newState = stateEst + kalmanGain * (measuredVal - temp(0));
            newStateCov = pEst - kalmanGain * M.transpose() * pEst;

            /*
            std::cout << "---------------beg" << std::endl;
            std::cout << "newState kalman:\n" << newState << std::endl;
            std::cout << "newStateCov:\n" << newStateCov << std::endl;
            std::cout << "state:\n" << state << std::endl;
            std::cout << "measuredVal:\n" << measuredVal << std::endl;
            std::cout << "stateCov:\n" << stateCov << std::endl;
            std::cout << "phi:\n" << phi << std::endl;
            std::cout << "Q:\n" << Q << std::endl;
            std::cout << "sigmaObs:\n" << sigmaObs << std::endl;
            std::cout << "M:\n" << M << std::endl;
            std::cout << "kalmanGain:\n" << kalmanGain << std::endl;
            std::cout << "stateEst:\n" << stateEst << std::endl;
            std::cout << "pEst:\n" << pEst << std::endl;
            std::cout << "-------------end" << std::endl;
            std::cin.get();
            */
            
        }

        float dataFs = 30000;
        const int fs;
    private:
        bool hasBeenUsed = false;
        bool freqsSet = false;

       // std::ofstream myfile;
        
        int histSize;
        int stride;

        // SSPE decs
        Matrix phi;
        Matrix Q;
        vector M;
       
        Array<Matrix> allXEB;

        Array<Matrix> allPEB;

        Array<float> dsBuf;
        Array<float> prevBuf;
        Array<Matrix> allX;
        Array<Matrix> allP;
        Array<std::complex<double>> complexArray;
        
        vector ampVec;
        vector sigmaFreqs;

        vector freqs;

        Matrix prevState;
        Matrix prevCov;

        // Initial Estimates
        vector ampEst;
        vector qEst;
        double obsErrorEst;


         
        int foiIndex;
        double sigmaObs;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SSPE);
    };
}

#endif 
